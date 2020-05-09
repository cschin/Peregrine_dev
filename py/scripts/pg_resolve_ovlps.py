import peregrine
from pypeflow.simple_pwatcher_bridge import (PypeProcWatcherWorkflow, Dist)
from pypeflow.tasks import gen_task
from docopt import docopt
import logging
import os
import sys

__version__ = peregrine.__version__

__doc__ = """
Workflow to resolve ovarlape due to reads are very similar to
each other either from different haplotypes or duplications.

Usage:
  pg_resolve_ovlps.py resolve <read_db_prefix> <ovlp_file>
                              <ovlp_ann_nchunk> <ovlp_ann_nproc>
                              <cns_nchunk> <cns_nproc>
                              --workdir <workdir>
                              [--min_len <min_len>]
                              [--min_idt <min_idt>]

  pg_resolve_ovlp.py (-h | --help)

  pg_resolve_ovlp.py --verison

Options:
  -h --help                   Show this help
  --version                   Show version
  --min_len <min_len>         Minimum overlap length for assembly graph construction [default: 4000]
  --min_idt <min_idt>         Minimum identity for considering two reads that are properly overlaps [default: 96]

"""

# Simple local-only submit-string.
submit = 'bash -C ${CMD} >| ${STDOUT_FILE} 2>| ${STDERR_FILE}'

def Task(script="", inputs={}, outputs={}, parameters=None, dist=None):
    if parameters is None:
        parameters = dict()
    if dist is None:
        dist = Dist()

    # Make paths relative to CWD. (But ok if caller does this.)
    def get_rel(maybe_abs):
        rel = dict()
        for (k, v) in maybe_abs.items():
            try:
                if os.path.isabs(v):
                    v = os.path.relpath(v)
                rel[k] = v
            except Exception:
                LOG.exception('Error for {!r}->{!r}'.format(k, v))
                raise
        return rel
    inputs = get_rel(inputs)
    outputs = get_rel(outputs)

    # All outputs must be in same directory.

    params = dict(parameters)

    pt = gen_task(script, inputs, outputs, params, dist)

    return pt

def run_overlap_annotation(wf, args):
    read_db_prefix = args["<read_db_prefix>"]
    ovlp_file = args["<ovlp_file>"]
    n_chunk = int(args["<ovlp_ann_nchunk>"])
    n_proc = int(args["<ovlp_ann_nproc>"])
    ovlp_ann_script = """
        /usr/bin/time -v shmr_overlap_annotate2 \
                -p {params.read_db_prefix}\
                -l {input.ovlp_file}\
                -t {params.n_chunk}\
                -c {params.my_chunk} > {params.ovlp_ann_fn} 
    """
    ovlp_dir = os.path.join(os.path.abspath(args["--workdir"]), "2-ovlp-rr")
    outputs = {}
    for my_chunk in range(1, n_chunk+1):
        ## each sub-job needs its own directory
        ovlp_ann_chunk_dir = os.path.join(ovlp_dir, f"chunk-{my_chunk:02d}")
        ovlp_ann_chunk_abs_prefix = os.path.join(ovlp_ann_chunk_dir, "ovlp.ann")
        ovlp_ann_fn = f"{ovlp_ann_chunk_abs_prefix}.{my_chunk:02d}"
        outputs[f'ovlp_ann_{my_chunk:02d}'] = ovlp_ann_fn
        
        wf.addTask(Task(
            script=ovlp_ann_script,
            inputs={
                'read_db': f"{read_db_prefix}.seqdb",
                'seqidx': f"{read_db_prefix}.idx",
                'ovlp_file': ovlp_file
            },
            outputs={f'ovlp_ann_{my_chunk:02d}': ovlp_ann_fn},
            parameters={
                'read_db_prefix': read_db_prefix,
                'n_chunk': n_chunk,
                'my_chunk': my_chunk,
                'ovlp_ann_fn': ovlp_ann_fn
            },
            dist=Dist(NPROC=1, local=True)
        ))

    wf.max_jobs = n_proc
    wf.refreshTargets()

    return outputs

def run_ovlp_to_ctg(wf, args, read_db_abs_prefix, read_db, ovlps):
    asm_script = """\

cat {params.ovlps} > preads-ann.ovl

cat preads-ann.ovl | awk '$NF == 1 {{print $1}}' | sort -u > filtered_out_reads
filter_ovlp.py > preads-rr.ovl
rm preads-ann.ovl
/usr/bin/time ovlp_to_graph.py  --overlap-file preads-rr.ovl --min_len {params.min_len} \
    --min_idt {params.min_idt} >& asm.log
/usr/bin/time graph_to_path.py >& to_path.log
/usr/bin/time path_to_contig.py {params.read_db_prefix} \
    p_ctg_tiling_path > {output.p_ctg} 2> to_contig.loga \


# TODO, change these parts using 'parallel' quikc hacks to use proper pypeflow tasks

mkdir -p cns_log
for i in `seq 1 {params.cns_nchunk}`; do 
   echo "/usr/bin/time tiling_path_ec.py {params.read_db_prefix} p_ctg_tiling_path preads-rr.ovl {params.cns_nchunk} $i  > p_cns_part_$i.fa 2> cns_log/log.p_cns.$i"; 
done | parallel -j {params.cns_nproc}
cat p_cns_part_*.fa > p_ctg_cns.fa
# rm p_cns_part_*.fa

for i in `seq 1 {params.cns_nchunk}`; do 
   echo "/usr/bin/time tiling_path_ec.py {params.read_db_prefix} a_ctg_tiling_path preads-rr.ovl {params.cns_nchunk} $i  > a_cns_part_$i.fa 2> cns_log/log.a_cns.$i"; 
done | parallel -j {params.cns_nproc}
cat a_cns_part_*.fa > a_ctg_cns.fa
# rm a_cns_part_*.fa
"""
    asm_dir = os.path.join(os.path.abspath(args["--workdir"]), "3-asm-rr")
    ovlps_list = " ".join(sorted([v for v in ovlps.values()]))
    inputs = {}
    inputs.update(ovlps)
    inputs.update(read_db)
    outputs = {}
    outputs["p_ctg"] = os.path.join(asm_dir, "p_ctg.fa")
    wf.addTask(Task(
        script=asm_script,
        inputs=inputs,
        outputs=outputs,
        parameters={
            'read_db_prefix': read_db_abs_prefix,
            'ovlps': ovlps_list,
            'min_len': int(args["--min_len"]),
            'min_idt': int(args["--min_idt"]),
            'cns_nchunk': int(args["<cns_nchunk>"]),
            'cns_nproc': int(args["<cns_nchunk>"])
        },
        dist=Dist(NPROC=1, local=True)
    ))
    wf.max_jobs = 1
    wf.refreshTargets()
    return outputs

def resolve_ovlps(args):

    job_defaults = dict(
        njobs=1,
        NPROC=1,
        MB=24000,
        submit=submit,
        job_type='local',
        pwatcher_type='blocking',
    )

    wf = PypeProcWatcherWorkflow(
            job_defaults=job_defaults,
    )
    print(args)
    ovlp_out = run_overlap_annotation(wf, args)
   
    read_db_abs_prefix = args["<read_db_prefix>"]
    
    read_db = {'read_db': f"{read_db_abs_prefix}.seqdb",
               'seqidx': f"{read_db_abs_prefix}.idx"}

    ctg_out = run_ovlp_to_ctg(wf, args,
                              read_db_abs_prefix,
                              read_db,
                              ovlp_out)



if __name__ == "__main__":
    import pkg_resources
    short_doc = """
Peregrine
=========

Peregrine is a fast genome assembler for accurate long
reads (length > 10kb, accuraccy > 99%). It can assemble
a human genome from 30x reads within 20 cpu hours from
reads to polished consensus. It uses Sparse HIereachical
MimiMizER (SHIMMER) for fast read-to-read overlaps without
explicitly quadratic comparisions used in other OLC
assemblers.

Peregrine Assembler and SHIMMER Genome Assembly Toolkit
Copyright (c) 2019- by Jason, Chen-Shan, Chin

Peregrine Assembler and  SHIMMER Genome Assembly Toolkit
is licensed under a Creative Commons
Attribution-NonCommercial-ShareAlike 4.0 International
License.

You should have received a copy of the license along with
this work. If not, see
<http://creativecommons.org/licenses/by-nc-sa/4.0/>.

************************************************************
If you want to use it for any commericial purposes
(including promotion activity for a commerical product),
please contact Jason Chin for a commericial license.
************************************************************

run `pg_resolve_ovlp.py -h` for help and other license information


"""
    sys.stderr.write(short_doc)


    logging.basicConfig(level=logging.INFO)
    args = docopt(__doc__, version=__version__)
    print(f"Peregrine Assembler ({__version__}) has been started with the following option:\n", args, file=sys.stderr)
    if args["resolve"]:
        resolve_ovlps(args)
