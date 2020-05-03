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

sig_pattern = re.compile(rb"[^-]{12}[-][^-]{12}")


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
    with open("preads.ovl") as f:
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
        hpc_seqs = {}
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
            hpc_seqs[r] = seq1

            aligend, aln, t_offset = get_align_range(seq1, seq0, 0,
                                                     max_dist=250,
                                                     aln_str=True)

            if aligend is not True:
                continue

            q_aln_str = falcon_ffi.string(aln.q_aln_str)
            t_aln_str = falcon_ffi.string(aln.t_aln_str)
            for m in sig_pattern.finditer(q_aln_str):
                s, e = m.span()
                p = q_aln_str[s:e], t_aln_str[s:e]
                sig_list.setdefault(p, [])
                sig_list[p].append(rid2)
            if aln is not None:
                falcon4py.free_alignment(aln)
            aln_set.add(rid2)

        sig_count = Counter(sig_list)

        marker_count = {}
        read_to_marker_count = {}
        for k, v in sig_count.items():
            m0, m1 = k
            m0_ = m0.replace(b"-", b"")
            m1_ = m1.replace(b"-", b"")
            c0 = 0
            c1 = 0
            for r in ovlps[r0]:
                read_to_marker_count.setdefault(r, [[], []])
                seq = hpc_seqs[r]
                tmp = len(re.findall(m0_, seq))
                if tmp > 0:
                    c0 += tmp
                    read_to_marker_count[r][0].append(m0)

                tmp = len(re.findall(m1_, seq))
                if tmp > 0:
                    c1 += tmp
                    read_to_marker_count[r][1].append(m1)
                    # print(m1)

            if c0 > 2 and c1 > 2:
                marker_count[m0] = c0
                marker_count[m1] = c1

        for r in ovlps[r0]:
            if r not in read_to_marker_count:
                c0 = 0
                c1 = 0
            else:
                c0 = 0
                for m in read_to_marker_count[r][0]:
                    if marker_count.get(m, 0) > 0:
                        c0 += 1
                    #print(0, m, marker_count.get(m, 0))
                #print()
                c1 = 0
                for m in read_to_marker_count[r][1]:
                    if marker_count.get(m, 0) > 0:
                        c1 += 1
                    #print(1, m, marker_count.get(m,0))
                #print(r, c0 ,c1)
                #print()
            flag = 0
            if c0 == 0 and c1 > 0:
                flag = 1
            d = ovlps[r0][r][1]

            print(" ".join(d), flag, c0, c1)


