from peregrine._falcon4py import ffi as falcon_ffi
from peregrine._falcon4py import lib as falcon4py
from peregrine.utils import SequenceDatabase
from peregrine.utils import rc
from peregrine.utils import get_align_range
from collections import Counter
import re, sys

## currently for interanl use
read_db_prefix = sys.argv[1]
total_chunks = int(sys.argv[2])
my_chunk = int(sys.argv[3])

sig_pattern = re.compile(rb"[^-]{12}-[^-]{12}")


def hp_compress(s):
    rtn = [s[0]]
    for c in s[1:]:
        if c != rtn[-1]:
            rtn.append(c)
    return bytes(rtn)



if __name__ == "__main__":

    ovlps = {}
    read_sdb = SequenceDatabase(
        f"{read_db_prefix}.idx",
        f"{read_db_prefix}.seqdb")
    with open("preads.ovl.2") as f:
        for r in f:
            r = r.strip().split()
            if r[0] == "-":
                continue
            id1, id2 = r[:2]
            if int(id1) % total_chunks == my_chunk:
                ovlps.setdefault(id1, {})
                ovlps[id1][id2] = (0, r)
            if int(id2) % total_chunks == my_chunk:
                ovlps.setdefault(id2, {})
                ovlps[id2][id1] = (1, r)

    for r0 in ovlps:
        if len(ovlps[r0]) < 5:
            continue
        rid = int(r0)
        read0 = read_sdb.get_subseq_by_rid(rid)

        sig_list = {}
        aln_set = set()
        for r in ovlps[r0]:
            #print(k,r, ovlps[k][r])
            rid = int(r0)
            d = ovlps[r0][r][1]
            if ovlps[r0][r][0] == 0:
                rid1 = int(d[0])
                rid2 = int(d[1])
                s1, e1, l1, direction, s2, e2, l2 = (int(_) for _ in d[5:12])
                e1 = e1 - 1
                e2 = e2 - 1
                assert e1 > s1
                if s1 < 0: s1 = 0
                if s2 < 0: s2 = 0
                seq0 = hp_compress(read0[s1:e1])
                seq1 = read_sdb.get_subseq_by_rid(rid2, s2, e2)
                if direction == 1:
                    seq1 = rc(seq1)
            else:
                rid2 = int(d[0])
                rid1 = int(d[1])
                s2, e2, l2, direction, s1, e1, l1 = (int(_) for _ in d[5:12])
                e1 = e1 - 1
                e2 = e2 - 1
                assert e1 > s1
                if s1 < 0: s1 = 0
                if s2 < 0: s2 = 0
                seq0 = hp_compress(read0[s1:e1])
                seq1 = read_sdb.get_subseq_by_rid(rid2, s2, e2)
                if direction == 1:
                    seq1 = rc(seq1)

            seq1 = hp_compress(seq1)

            aligend, aln, t_offset = get_align_range(seq1, seq0, 0,
                                                     max_dist=250,
                                                     aln_str=True)

            if aligend is not True:
                continue

            s = falcon_ffi.string(aln.q_aln_str)
            for p in sig_pattern.findall(s):
                sig_list.setdefault(p, [])
                sig_list[p].append(rid2)
            if aln is not None:
                falcon4py.free_alignment(aln)
            aln_set.add(rid2)

        sig_count = Counter(sig_list)
        anti_phase_set = set()
        for k, v in sig_count.items():
            if abs(len(v)/len(aln_set) - 0.5) > 0.4:
                continue
            if len(v) < 3 or len(aln_set)-len(v) < 3:
                continue
            if r0 not in set(v):
                for rid in v:
                    assert r0 != rid
                    anti_phase_set.add(rid)

        #anti_phase_set, aln_set = get_aln_signatures(seqs, rid)
        #print(anti_phase_set, aln_set, len(anti_phase_set), len(aln_set), len(ovlps[r0]))
        for r in ovlps[r0]:
            if int(r) in anti_phase_set:
                flag = 0
            else:
                flag = 1
            d = ovlps[r0][r][1]
            print(" ".join(d), flag)

