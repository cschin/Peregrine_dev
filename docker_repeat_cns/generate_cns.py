
from collections import Counter
from collections import defaultdict
import random
import glob
import sys
import os

import peregrine.utils
from peregrine.utils import SequenceDatabase
from peregrine.utils import get_shimmers_from_seq
from peregrine.utils import get_shimmer_alns_from_seqs
from peregrine.utils import rc
from peregrine.utils import get_cns_from_reads

import umap
import numpy as np
from sklearn.cluster import DBSCAN





def hp_compress(s):
    rtn = [s[0]]
    pos = 0
    coor_map = [pos]
    for c in s[1:]:
        if c != rtn[-1]:
            rtn.append(c)
            pos += 1
        coor_map.append(pos)
    assert len(coor_map) == len(s)
    return bytes(rtn), coor_map

def get_signature_mers(compressed_reads, kmer_size = 20):
    signature_mers = defaultdict(int)

    for rid in compressed_reads:
        c = defaultdict(int)
        s = compressed_reads[rid]
        kmers = get_shimmers_from_seq(s, levels=1, reduction_factor=1, w=kmer_size, k=kmer_size)
        kmers = [kmers[j].mmer for j in range(len(kmers))]
        for mer in kmers:
            c[mer] += 1

        for mer in kmers:
            if c[mer] == 1:
                signature_mers[mer] += 1
    signature_mers2 = {}
    for mer in signature_mers:
        c = signature_mers[mer]
        if c < len(compressed_reads)*0.2 or c > len(compressed_reads)*0.8:
            continue
        signature_mers2[mer] = c
    return signature_mers, signature_mers2


def get_signature_vectors(compressed_reads, signature_mers, kmer_size = 20):
    read2path = {}
    for rid in compressed_reads:
        c = defaultdict(int)
        s = compressed_reads[rid]
        path = []
        kmers = get_shimmers_from_seq(s, levels=1, reduction_factor=1, w=kmer_size, k=kmer_size)
        kmers = [kmers[j].mmer for j in range(len(kmers))]
        for mer in kmers:
            if mer in signature_mers:
                path.append(mer)
        read2path[rid] = path
    index=dict(zip(signature_mers.keys(), range(len(signature_mers))))
    v = {}
    for i in read2path:
        vv = np.zeros(len(index))
        for m in set(read2path[i]):
            vv[index[m]]=1

        v[i] = vv
    v2 = np.array(list(v.values()))
    return v2


def get_consensus(rids, labels, read_sdb):
    consensus = []
    for cl in Counter(labels):
        if cl == -1:
            continue
        reads = {}
        read_length = {}
        for rid in rids[labels==cl]:
            s = read_sdb.get_subseq_by_rid(rid)
            reads[rid] = s
            read_length[rid] = len(s)
        #print(len(reads), np.mean(read_length), np.max(read_length), np.min(read_length), np.median(read_length))
        reads2 = {}
        median = np.median(list(read_length.values()))
        #c = 0
        for rid, s in reads.items():
            if abs(read_length[rid] - median) < 500:
                reads2[rid] = s
                #c += 1
        #plt.figure()
        #plt.hist(list(read_length.values()), bins=32, range=(10000,15000))
        cns=get_cns_from_reads(list(reads2.values()), sort_reads=False, best_n=1000, levels=1, w=32)
        consensus.append( (len(reads2), cns, reads, reads2))
    consensus.sort()
    consensus.reverse()
    return consensus


def prepare_reads(workdir, bamfile, ref_file, prefix=None):
    if prefix is None:
        prefix = os.path.basename(bamfile)
        prefix = ".".join(prefix.split(".")[:-2])

    os.system(f"mkdir -p {workdir}")
    os.system(f"samtools fastq  -0 {workdir}/{prefix}.fq {bamfile}")
    os.system(f"echo {workdir}/{prefix}.fq > {workdir}/reads.lst")
    os.system(f"shmr_mkseqdb -d {workdir}/reads.lst -p {workdir}/reads")
    os.system(f"echo {ref_file} > {workdir}/ref.lst")
    os.system(f"shmr_mkseqdb -d {workdir}/ref.lst -p {workdir}/ref")


if __name__ == "__main__":
    kmer_size = 20
    random.seed(2501)

    workdir = sys.argv[1]
    outdir = sys.argv[2]
    bamfile = sys.argv[3]
    reffile = sys.argv[4]
    min_len = int(sys.argv[5])
    max_len = int(sys.argv[6])

    os.system(f"mkdir -p {outdir}")
    prepare_reads(workdir, bamfile, reffile)

    prefix = os.path.basename(bamfile)
    prefix = ".".join(prefix.split(".")[:-2])

    read_sdb=SequenceDatabase(f"{workdir}/reads.idx", f"{workdir}/reads.seqdb")
    sub_ref_sdb=SequenceDatabase(f"{workdir}/ref.idx", f"{workdir}/ref.seqdb")

    rids = list([rid for rid in read_sdb.index_data.keys()
                 if read_sdb.index_data[rid].length > min_len and
                    read_sdb.index_data[rid].length < max_len])
    if len(rids) > 500:
        rids = np.array(sorted(random.sample(rids, 500)))

    compressed_reads= {}
    for rid in rids:
        seq = read_sdb.get_subseq_by_rid(rid)
        s, _ = hp_compress(seq)
        compressed_reads[rid] = s

    rids = np.array(rids)
    print(len(compressed_reads))
    mer, mer2 = get_signature_mers(compressed_reads, kmer_size = kmer_size)
    v = get_signature_vectors(compressed_reads, mer2, kmer_size = kmer_size)
    print(v.shape)
    embedding = umap.UMAP(n_neighbors=20,
                      min_dist=0.5,
                      random_state=42,
                      metric='hamming').fit_transform(v)

    clustering = DBSCAN(eps=1, min_samples=20).fit(embedding)

    labels = clustering.labels_
    consensus = get_consensus(rids, labels, read_sdb)

    s = sub_ref_sdb.get_subseq_by_rid(0)

    cns_stats = {}
    i = 1
    summary_file = open(f"{outdir}/{prefix}.summary","w")

    for c, s2, reads, reads2 in consensus:
        r = c / len(compressed_reads) * 100
        m = len(compressed_reads)
        len_ = len(s2)
        q = len(reads2)/len(reads)*100.0
        fnprefix = f"{outdir}/{prefix}.cns{i:02d}"
        fn = f"{fnprefix}.fa"
        with open(fn, "w") as f:
            print(f">{prefix}-cns{i:02d}  n={c} m={m} r={r:0.1f} q={q:0.1f} l={len_}", file=summary_file)
            print(f">{prefix}-cns{i:02d}  n={c} m={m} r={r:0.1f} q={q:0.1f} l={len_}", file=f)
            print(s2.decode(), file=f)

        cmds = [f"/opt/dipcall.kit/minimap2 -a -xasm5 --cs -r2k -z1000000,100000 -t8 {reffile} {fn} 2> {fnprefix}.sam.gz.log | gzip > {fnprefix}.sam.gz"]
        cmds.append(f"/opt/dipcall.kit/k8 /opt/dipcall.kit/dipcall-aux.js samflt -L 10000 {fnprefix}.sam.gz | /opt/dipcall.kit/samtools sort -m4G --threads 1 -o {fnprefix}.bam -")
        cmds.append(f"/opt/dipcall.kit/htsbox pileup -q5 -evcf {reffile} {fnprefix}.bam  | /opt/dipcall.kit/htsbox bgzip > {fnprefix}.vcf.gz")
        cmds.append(f"rm {fnprefix}.sam.gz.log {fnprefix}.sam.gz {fnprefix}.bam ")
        print("\n".join(cmds))
        print()
        os.system("\n".join(cmds))

        cns_stats[f"{prefix}.cns{i:02d}"] = (c, r, len_)
        fn = f"{outdir}/reads-{prefix}.cns{i:02d}.fa"
        with open(fn, "w") as f2:
            for rid in reads2:
                print(f">{rid}", file=f2)
                print(reads[rid].decode(), file=f2)

        i+=1
    summary_file.close()


