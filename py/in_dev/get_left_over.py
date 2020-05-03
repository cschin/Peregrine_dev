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


ctg_sdb=SequenceDatabase("pa_phase_ctg.idx", "pa_phase_ctg.seqdb")

htig_paths = {}
ctgs = set()
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
        if source != "backbone":
            ctgs.add(source)
        ss = int(r[6])
        st = int(r[7])
        direction = int(r[8])
        htig_paths.setdefault(htig, [])
        htig_paths[htig].append( (backbone, bs, bt, source, ss, st, direction) )


for rid in ctg_sdb.index_data:
    ctg_name = ctg_sdb.index_data[rid].rname
    if ctg_name not in ctgs:
        print(f">{ctg_name}")
        seq = ctg_sdb.get_subseq_by_rid(rid)
        print(seq.decode())
