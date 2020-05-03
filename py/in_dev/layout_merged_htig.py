import peregrine.utils
from peregrine.utils import SequenceDatabase
from peregrine.utils import get_shimmers_from_seq
from peregrine.utils import get_shimmer_alns
from peregrine.utils import get_cigar
from peregrine.utils import get_shimmer_alns_from_seqs
from peregrine.utils import rc
from peregrine.ctg_align import SeqDBAligner, get_shimmer_dots
from intervaltree import Interval, IntervalTree
from collections import Counter
import networkx as nx
import numpy as np
from multiprocessing import Pool


ref_sdb=SequenceDatabase("p_ctg_cns.idx", "p_ctg_cns.seqdb")
ctg_sdb=SequenceDatabase("pa_phase_ctg.idx", "pa_phase_ctg.seqdb")

htig_paths = {}
with open("tmp") as f:
    for r in f:
        r = r.strip().split()
        """
        000001F-H1 000001F 1077129 1080126 2997 backbone 1077129 1080126 0
        """
        htig = r[0]
        backbone = r[1]
        bs = int(r[2])
        bt = int(r[3])
        source = r[5]
        ss = int(r[6])
        st = int(r[7])
        direction = int(r[8])
        htig_paths.setdefault(htig, [])
        htig_paths[htig].append( (backbone, bs, bt, source, ss, st, direction) )


def get_new_seq(input_):
    global ref_sdb, ctg_sdb
    htig, htig_paths = input_ 
    seqs = []
    backbone, bs, bt, source, ss, st, direction = htig_paths[htig][0]
    if source == "backbone":
        seq = ref_sdb.get_subseq_by_name(backbone, 0, bt)
    else:
        seq = ctg_sdb.get_subseq_by_name(source, 0, st, direction)
    seqs.append(seq)
    for seg in htig_paths[htig][1:-1]:
        backbone, bs, bt, source, ss, st, direction = seg
        if source == "backbone":
            seq = ref_sdb.get_subseq_by_name(backbone, bs, bt)
        else:
            seq = ctg_sdb.get_subseq_by_name(source, ss, st, direction)
        seqs.append(seq)
    if len(htig_paths[htig]) >= 2:
        backbone, bs, bt, source, ss, st, direction = htig_paths[htig][-1]
        if source == "backbone":
            slen = ref_sdb.get_seq_index_by_name(backbone).length
            seq = ref_sdb.get_subseq_by_name(backbone, bs, slen )
        else:
            slen = ctg_sdb.get_seq_index_by_name(source).length
            seq = ctg_sdb.get_subseq_by_name(source, ss, slen, direction)
        seqs.append(seq)
    out_seq = b"".join(seqs).decode()
    return htig, out_seq

htig_names = sorted(htig_paths.keys())
h1file = open("H1.fa", "w")
h2file = open("H2.fa", "w")

input_list = []
for htig in htig_names:
    input_list.append( (htig, htig_paths) ) 

p = Pool(20)
for htig, seq in p.imap(get_new_seq, input_list):
    if htig.split("-")[-1] == "H1":
        outfile = h1file
    else:
        outfile = h2file

    print(f">{htig}", file=outfile)
    print(seq, file=outfile)

h1file.close()
h2file.close()

