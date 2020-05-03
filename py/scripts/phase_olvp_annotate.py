from peregrine._falcon4py import ffi as falcon_ffi
from peregrine._falcon4py import lib as falcon4py
from peregrine.utils import SequenceDatabase
from peregrine.utils import rc
from peregrine.utils import get_align_range
from collections import defaultdict
import numpy as np
import sys

## currently for interanl use
read_db_prefix = sys.argv[1]
overlap_file = sys.argv[2]
total_chunks = int(sys.argv[3])
my_chunk = int(sys.argv[4])



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


def get_kmer_count(rid, ovlps, k=31):
    rseq = read_sdb.get_subseq_by_rid(rid)
    rseq, coor_map = hp_compress(rseq)
    cov = np.zeros(len(rseq)-k)
    count = defaultdict(int)
    support_reads = {}
    for rid2 in ovlps[rid]:
        #print(ovlps[rid][rid2])
        rec = ovlps[rid][rid2][0]
        direction = int(ovlps[rid][rid2][0][8])
        if ovlps[rid][rid2][1] == 0:
            ref_s = int(rec[5])
            ref_e = int(rec[6])
            s = int(rec[9])
            e = int(rec[10])
            l = int(rec[11])
            if direction == 1:
                s, e = l - e, l - s
            seq = read_sdb.get_subseq_by_rid(rid2, s, e, direction=direction)
        else:
            ref_s = int(rec[9])
            ref_e = int(rec[10])
            ref_l = int(rec[11])
            s = int(rec[5])
            e = int(rec[6])
            l = int(rec[7])
            if direction == 1:
                s, e = l - e, l - s
            seq = read_sdb.get_subseq_by_rid(rid2, s, e, direction=direction)
        
        seq, _ = hp_compress(seq)
        support_reads[rid2] = seq, coor_map[ref_s], coor_map[ref_e-1], direction, ovlps[rid][rid2][1]
        cov[coor_map[ref_s]:coor_map[ref_e-1]-k+1] += 1 
        for i in range(0, len(seq)-k):
            count[seq[i:i+k]] += 1

    c = np.zeros(len(rseq)-k)
    for i in range(0, len(rseq)-k):
        c[i] = count.get(rseq[i:i+k], 0)

    markers = []
    if any(cov == 0):
        return None, None, None, None

    for i, r in enumerate(c/(cov+1)):
        if c[i] > 3 and r < 0.8:
            markers.append( (i, rseq[i:i+k], c[i], cov[i], r ))
            
    return count, support_reads, cov, markers

if __name__ == "__main__":

    ovlps = {}
    read_sdb = SequenceDatabase(
        f"{read_db_prefix}.idx",
        f"{read_db_prefix}.seqdb")

    with open(overlap_file) as f:
        for r in f:
            r = r.strip().split()
            if r[0] == "-":
                continue

            id1, id2 = r[:2]
            id1 = int(id1)
            id2 = int(id2)
            if id1 % total_chunks == my_chunk: 
                ovlps.setdefault(id1, {})
                ovlps[id1][id2] = (r, 0)
            if id2 % total_chunks == my_chunk: 
                ovlps.setdefault(id2, {})
                ovlps[id2][id1] = (r, 1)
    k = 31
    for rid in ovlps:
        count, support_reads, cov, markers = get_kmer_count(rid, ovlps, k=k)
        if count is None:
            continue
        for rid1 in ovlps[rid]:
            if ovlps[rid][rid1][1] == 1:
                continue
            seq1, ref_s, ref_e, direction, d = support_reads[rid1]
           
            m = set()
            for i, kmer, c0, c1, r in markers:
                if i >= ref_s and i < ref_e:
                    m.add(kmer)
            m2 = set()
            for i in range(len(seq1)-k):
                if seq1[i:i+k] in m:
                    m2.add(seq1[i:i+k])
            if len(m) == 0:
                print(" ".join(ovlps[rid][rid1][0]), len(m), len(m2), len(m)-len(m2), -1)
                continue
            print(" ".join(ovlps[rid][rid1][0]), len(m), len(m2), len(m)-len(m2), len(m2)/len(m)) 


