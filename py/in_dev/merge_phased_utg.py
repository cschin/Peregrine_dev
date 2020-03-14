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



def get_ctg2match(map_fn, ctg_sdb):
    ctg2match = {}
    rctg2path = {}
    with open(map_fn) as f:
        for row in f:
            row = [int(c) for c in row.strip().split()]
            ctg_name = ctg_sdb.index_data[row[3]].rname
            ctg2match.setdefault(ctg_name,[])
            ctg2match[ctg_name].append(row)
            
            rctg_name = ref_sdb.index_data[row[0]].rname 
            rctg2path.setdefault(rctg_name, {})
            k = (row[1], row[2])
            rctg2path[rctg_name].setdefault(k,[])
            rctg2path[rctg_name][k].append(row) 
    return ctg2match, rctg2path


def get_best_match(ctg2match):
    ctg2best = {}
    r_ctg2ctg = {}
    for ctg in ctg2match:
        c=Counter([(_[0], _[6]) for _ in ctg2match[ctg]])
        r_ctg_id, direction = c.most_common(1)[0][0]
        r_ctg = ref_sdb.index_data[r_ctg_id].rname
        ctg2best[ctg] = r_ctg, direction
        r_ctg2ctg.setdefault(r_ctg, [])
        r_ctg2ctg[r_ctg].append((ctg, direction))
    return ctg2best, r_ctg2ctg


def get_match_group(matches):
    groups = []
    group = []
    pos = -1
    pos2 = -1
    for m in matches:
        if pos == -1:
            pos = m[1]
            pos2 = m[4]
            group.append(m)
            continue
        c_pos = m[1]
        c_pos2 = m[4]
        if c_pos - pos < 0 or c_pos - pos > 50000:
            # print(c_pos, pos)
            groups.append(group)
            group = []
        if c_pos2 > pos2:
            group.append(m)
        pos = c_pos
        pos2 = c_pos2
    if len(group) != 0:
        groups.append(group)
    return groups

def get_ctg2group(r_ctg2ctg, rctg):
    ctg2groups = {}
    for ctg, direction in r_ctg2ctg[rctg]:
        #print(ctg)
        r_ctg_id = ref_sdb.name2rid[rctg]
        matches = []
        for r in ctg2match[ctg]:
            if r[0] != r_ctg_id or r[6] != direction:
                continue
            matches.append(r)
        m_group = get_match_group(matches)
        for mg in m_group:
            if len(mg) < 2:
                continue
            # print(ctg, matches[0], matches[-1], sep="\n")
            s0 = mg[0][1]
            t0 = mg[-1][2]
            s1 = mg[0][4]
            t1 = mg[-1][5]
            d1 = mg[0][6]
            d2 = mg[-1][6]
            ctg2groups.setdefault(ctg, [])
            ctg2groups[ctg].append( (s0, t0, s1, t1, direction) )
    return ctg2groups

def get_layout_segements(r_ctg2ctg, rctg, ctg_sdb, ctg2groups, phase="H1", exclude_ctgs = set()):
    used_ctgs = set() 
    segs = []
    seg_end = -1
    bb_bgn = -1
    bb_end = -1 
    for k in sorted(rctg2path[rctg].keys()):
        if k[0] < seg_end:
            continue
        
        candidates = []

        for m in rctg2path[rctg][k]:
            ctg = ctg_sdb.index_data[m[3]].rname
            if ctg not in ctg2groups:
                continue
            if ctg in exclude_ctgs:
                continue
            for g in ctg2groups[ctg]:
                candidates.append((g, ctg))
        
        if len(candidates) != 0:
            candidates2 = []
            for c in candidates:
                overlap = c[0][0] - seg_end
                extend = c[0][1] - seg_end
                if overlap > 10000:
                    continue
                    
                # if extend > 2*abs(overlap) and overlap > -50000:
                if extend > 2*abs(overlap):
                    candidates2.append(c)

            if len(candidates2) > 0:
                candidates2.sort(key=lambda x: x[0][1])
                c = candidates2[-1]
                overlap = c[0][0] - seg_end
                extend = c[0][1] - seg_end
                s0, t0, s1, t1, direction = c[0]
                ctg = c[1]
                for m in rctg2path[rctg][k]:
                    if ctg_sdb.index_data[m[3]].rname == c[1] and \
                       m[1] == k[0]:
                        s0 = m[1]
                        s1 = m[4]

                if bb_bgn == -1:
                    bb_bgn = t0
                else:
                    segs.append( (f"{rctg}-{phase}", f"{rctg}", bb_bgn, s0, s0-bb_bgn, "backbone", bb_bgn, s0, 0) )
                    bb_bgn = t0
                segs.append( (f"{rctg}-{phase}", f"{rctg}", s0, t0, t0-s0, f"{ctg}", s1, t1, direction) )
                used_ctgs.add(ctg)
                seg_end = c[0][1]
                
            else:
                if seg_end != k[1]:
                    seg_end = k[1]
    return segs, used_ctgs


ref_sdb=SequenceDatabase("p_ctg_cns.idx", "p_ctg_cns.seqdb")
ctg_sdb=SequenceDatabase("pa_phase_ctg.idx", "pa_phase_ctg.seqdb")

ctg2match, rctg2path = get_ctg2match("p_ctg_cns-pa_phase_ctg-map", ctg_sdb)
ctg2best, r_ctg2ctg = get_best_match(ctg2match)

all_used_ctgs = set()
f = open("ctg_assign.txt", "w")
for rctg in rctg2path:

    if rctg not in r_ctg2ctg:
        s0 = 0
        t0 = ref_sdb.get_seq_index_by_name(rctg).length
        segs = ((f"{rctg}-H1", f"{rctg}", s0, t0, t0-s0, "backbone", s0, t0, 0),)
    else:
        ctg2groups = get_ctg2group(r_ctg2ctg, rctg)
        segs, p_utgs = get_layout_segements(r_ctg2ctg, rctg, ctg_sdb,
                                               ctg2groups, phase="H1") 

    for seg in segs:
        print( " ".join([str(_) for _ in seg]) )

    if rctg not in r_ctg2ctg:
        s0 = 0
        t0 = ref_sdb.get_seq_index_by_name(rctg).length
        segs = ((f"{rctg}-H2", f"{rctg}", s0, t0, t0-s0, "backbone", s0, t0, 0),)
    else:
        ctg2groups = get_ctg2group(r_ctg2ctg, rctg)
        segs, a_utgs = get_layout_segements(r_ctg2ctg, rctg, ctg_sdb,
                                               ctg2groups, phase="H2", exclude_ctgs=p_utgs) 

    for seg in segs:
        print( " ".join([str(_) for _ in seg]) )

    for n in p_utgs:
        print(rctg, n, 0, file=f)

    for n in a_utgs:
        print(rctg, n, 1, file=f)
    
    if rctg in r_ctg2ctg:
        s = set([_[0] for _ in r_ctg2ctg[rctg]])
        for n in s - p_utgs - a_utgs:
            print(rctg, n, 2, file=f)


f.close()
