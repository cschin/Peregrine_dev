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

p_ctg_file = open("p_utg_test.fa", "w")
a_ctg_file = open("a_utg_test.fa", "w")
extra_ctg_file = open("extra_utg_test.fa", "w")
with open("ctg_assign.txt") as f:
    for r in f:
        r = r.strip().split()
        p_ctg, a_ctg, t = r
        seq = ctg_sdb.get_subseq_by_name(a_ctg)
        if int(t) == 0:
            print(f">{p_ctg}/{a_ctg}", file=p_ctg_file)
            print(seq.decode(), file=p_ctg_file)
        if int(t) == 1:
            print(f">{p_ctg}/{a_ctg}", file=a_ctg_file)
            print(seq.decode(), file=a_ctg_file)
        elif int(t) == 2:
            if len(seq) < 50000:
                continue
            print(f">{p_ctg}/{a_ctg}", file=extra_ctg_file)
            print(seq.decode(), file=extra_ctg_file)
p_ctg_file.close()
a_ctg_file.close()
extra_ctg_file.close() 
