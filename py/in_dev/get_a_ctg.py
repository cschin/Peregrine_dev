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
with open("ctg_assign.txt") as f:
    for r in f:
        r = r.strip().split()
        p_ctg, a_ctg, t = r
        if int(t) == 0:
            continue
        seq = ctg_sdb.get_subseq_by_name(a_ctg)
        if len(seq) < 50000:
            continue
        print(f">{p_ctg}/{a_ctg}")
        print(seq.decode())


