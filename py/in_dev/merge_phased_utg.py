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
    with open(map_fn) as f:
        for row in f:
            row = [int(c) for c in row.strip().split()]
            ctg_name = ctg_sdb.index_data[row[3]].rname
            ctg2match.setdefault(ctg_name,[])
            ctg2match[ctg_name].append(row)
    return ctg2match


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


def get_marker_link(rctg, ref_sdb, ctg2match, r_ctg2ctg, excluded=set()):
    markers = set()
    links = {}
    s0=ref_sdb.get_subseq_by_name(rctg)
    if rctg not in r_ctg2ctg:
        markers = (0, ref_sdb.get_seq_index_by_name(rctg).length-32)
        return markers, links

    for ctg, direction in r_ctg2ctg[rctg]:
        if ctg in excluded:
            continue
        r_ctg_id = ref_sdb.name2rid[rctg]
        matches = []
        for r in ctg2match[ctg]:
            if r[0] != r_ctg_id or r[6] != direction:
                continue
            #if r[-1] > 10:
            #    continue
            matches.append(r)
        # print(ctg, direction, len(matches))
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
            assert d1 == direction
            assert d2 == direction
            links[(s0, t0)] = ctg, s1, t1, direction
            markers.add(s0)
            markers.add(t0)
    if len(markers) < 2:
        links = {}
        markers = (0, ref_sdb.get_seq_index_by_name(rctg).length-32)
        return markers, links

    markers = list(markers)
    markers.sort()
    return markers, links

def get_filtered_links(rctg, links):
    intervals = []
    filtered_links = {}
    filtered_markers = set()
    
    keys = list(links.keys())
    keys.sort(key=lambda x:x[0]-x[1])

    for s, t in keys:
        flag = 0
        for itvl in intervals:
            if s > itvl[0] and s < itvl[1]:
                flag = 1
                continue
            if t > itvl[0] and t < itvl[1]:
                flag = 1
                continue
        if flag == 1:
            continue
        filtered_links[(s,t)] = links[(s,t)]
        intervals.append((s, t))
        filtered_markers.add(s)
        filtered_markers.add(t)
    filtered_markers = list(filtered_markers)
    filtered_markers.sort()
    if len(filtered_markers) < 2:
        links = {}
        markers = (0, ref_sdb.get_seq_index_by_name(rctg).length-32)
        return markers, links
    return filtered_markers, filtered_links


ref_sdb=SequenceDatabase("p_ctg_cns.idx", "p_ctg_cns.seqdb")
ctg_sdb=SequenceDatabase("pa_phase_ctg.idx", "pa_phase_ctg.seqdb")

ctg2match = get_ctg2match("p_ctg_cns-pa_phase_ctg-map", ctg_sdb)
ctg2best, r_ctg2ctg = get_best_match(ctg2match)


def get_paths(rctg_id):
    rtn = []
    rctg_name = ref_sdb.index_data[rctg_id].rname
    markers, links = get_marker_link(rctg_name, ref_sdb, ctg2match, r_ctg2ctg)
    markers, links = get_filtered_links(rctg_name, links)

    G=nx.DiGraph()
    nx.add_path(G, markers, ctg="backbone")
    for s, t in links:
        ctg=links[(s,t)][0]
        G.add_edge(s, t, ctg=ctg)
    #print(len(G))
    #nx.write_gexf(G,"test2.gexf")

    #if len(markers) < 2:
    #    print("XX", rctg_id, markers)
    path=nx.shortest_path(G, markers[0], markers[-1])
    edges = zip(path[0:-1], path[1:])
    main_ctgs = set()
    for s, t in edges:
        if (s, t) in links:
            ctg, s2, t2, d = links[(s,t)]
            rtn.append((f"{rctg_name}-H1", rctg_name, s, t, t-s, ctg, s2, t2, d))
            main_ctgs.add(ctg)
        else:
            rtn.append((f"{rctg_name}-H1", rctg_name, s, t, t-s, G[s][t]["ctg"], s, t, 0))

    markers, links = get_marker_link(rctg_name, ref_sdb, ctg2match, r_ctg2ctg, excluded=main_ctgs)
    markers, links = get_filtered_links(rctg_name, links)

    G=nx.DiGraph()
    nx.add_path(G, markers, ctg="backbone")
    for s, t in links:
        ctg=links[(s,t)][0]
        G.add_edge(s, t, ctg=ctg)
    #print(len(G))
    #nx.write_gexf(G,"test2.gexf")

    #if len(markers) < 2:
    #    print("XX", rctg_id, markers)
    path=nx.shortest_path(G, markers[0], markers[-1])
    edges = zip(path[0:-1], path[1:])
    secondary_ctgs = set()
    for s, t in edges:
        if (s, t) in links:
            ctg, s2, t2, d = links[(s,t)]
            rtn.append((f"{rctg_name}-H2", rctg_name, s, t, t-s, ctg, s2, t2, d))
            secondary_ctgs.add(ctg)
        else:
            rtn.append((f"{rctg_name}-H2", rctg_name, s, t, t-s, G[s][t]["ctg"], s, t, 0))
    additional_ctgs = set() 
    if rctg_name  in r_ctg2ctg:
        for ctg, _ in r_ctg2ctg[rctg_name]:
            if ctg not in main_ctgs and ctg not in secondary_ctgs:
                additional_ctgs.add(ctg)
    return rtn, rctg_name, main_ctgs, secondary_ctgs, additional_ctgs

p = Pool(20)
f = open("ctg_assign.txt", "w")
for rtn, rctg_name, main_ctgs, secondary_ctgs, additional_ctgs in p.imap(get_paths, ref_sdb.index_data.keys()):
    for l in rtn:
        print(" ".join([str(_) for _ in l]))
    all_ctgs = main_ctgs | secondary_ctgs | additional_ctgs
    for ctg in list(all_ctgs):
        if ctg in main_ctgs:
            print( rctg_name, ctg, 0, file = f)
        elif ctg in secondary_ctgs:
            print( rctg_name, ctg, 1, file = f)
        elif ctg in additional_ctgs: 
            print( rctg_name, ctg, 2, file = f)
f.close()

#for rtn in map(get_paths, ref_sdb.index_data.keys()):
#    for l in rtn:
#        print(" ".join([str(_) for _ in l]))

