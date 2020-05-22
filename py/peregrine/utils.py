# import sys
# import os
import numpy as np
import mmap
import sys
from collections import namedtuple
from ._shimmer4py import ffi as shimmer_ffi
from ._shimmer4py import lib as shimmer4py
from ._falcon4py import ffi as falcon_ffi
from ._falcon4py import lib as falcon4py
from collections import Counter
# from ._ksw4py import ffi as ksw4py_ffi


rmap = dict(list(zip(b"ACGT", b"TGCA")))


def rc(seq):
    return bytes([rmap[c] for c in seq[::-1]])


MMER = namedtuple("MMER", "mmer span rid, pos_end, direction")


def mmer2tuple(mmer):
    x = mmer.x
    y = mmer.y
    span = x & 0xFF
    mmer = x >> 8
    rid = y >> 32
    pos_end = ((y & 0xFFFFFFFF) >> 1) + 1
    direction = y & 0x1
    return MMER(mmer, span, rid, pos_end, direction)


class Shimmer(object):
    def __init__(self):
        self.mmers = shimmer_ffi.new("mm128_v *")

    def __len__(self):
        return self.mmers.n

    def __getitem__(self, i):
        assert i < self.mmers.n
        return mmer2tuple(self.mmers.a[i])

    def __del__(self):
        shimmer4py.free(self.mmers.a)
        shimmer_ffi.release(self.mmers)


def get_shimmers_from_seq(seq, rid=0,
                          levels=2, reduction_factor=3,
                          k=16, w=80):
    assert levels <= 2
    c_null = shimmer_ffi.NULL
    mmers_L0 = Shimmer()
    shimmer4py.mm_sketch(c_null, seq, len(seq), w, k, rid, 0, mmers_L0.mmers)
    if levels == 0:
        return mmers_L0
    elif levels == 1:
        mmers_L1 = Shimmer()
        shimmer4py.mm_reduce(mmers_L0.mmers, mmers_L1.mmers, reduction_factor)
        del mmers_L0
        return mmers_L1
    elif levels == 2:
        mmers_L1 = Shimmer()
        mmers_L2 = Shimmer()
        shimmer4py.mm_reduce(mmers_L0.mmers, mmers_L1.mmers, reduction_factor)
        shimmer4py.mm_reduce(mmers_L1.mmers, mmers_L2.mmers, reduction_factor)
        del mmers_L1
        del mmers_L0
        return mmers_L2


ALIGNED_MMER = namedtuple("ALIGNED_MMER", "mmer0 mmer1")


def get_shimmer_alns(shimmers0, shimmers1,
                     direction=0,
                     max_diff=100,
                     max_dist=1200,
                     max_repeat=1):
    _shimmers0 = shimmers0.mmers
    _shimmers1 = shimmers1.mmers
    aln = shimmer4py.shmr_aln(_shimmers0, _shimmers1, direction,
                              max_diff, max_dist, max_repeat)
    aln_chains = []
    for i in range(aln.n):
        chain = []
        offsets = np.zeros(aln.a[i].idx0.n, dtype=np.float)
        for j in range(aln.a[i].idx0.n):
            idx0, idx1 = aln.a[i].idx0.a[j], aln.a[i].idx1.a[j]
            mmer0 = mmer2tuple(_shimmers0.a[idx0])
            mmer1 = mmer2tuple(_shimmers1.a[idx1])
            chain.append(ALIGNED_MMER(mmer0, mmer1))
            if direction == 0:  # same direction
                d = mmer0[3] - mmer1[3]
            else:
                d = mmer0[3] + mmer1[3]
            offsets[j] = d
        aln_chains.append((chain, (np.max(d), np.mean(d), np.min(d))))
    shimmer4py.free_shmr_alns(aln)
    return aln_chains


def get_shimmer_alns_from_seqs(seq0, seq1, parameters={}):
    reduction_factor = parameters.get("reduction_factor", 3)
    direction = parameters.get("direction", 0)
    max_diff = parameters.get("max_diff", 1000)
    max_dist = parameters.get("max_dist", 15000)
    max_repeat = parameters.get("max_repeat", 1)
    w = parameters.get("w", 80)
    k = parameters.get("k", 16)

    seq0_shimmers = get_shimmers_from_seq(
        seq0,
        rid=0,
        reduction_factor=reduction_factor,
        k=k, w=w)

    seq1_shimmers = get_shimmers_from_seq(
        seq1,
        rid=1,
        reduction_factor=reduction_factor,
        k=k, w=w)

    shimmer_alns = get_shimmer_alns(seq0_shimmers,
                                    seq1_shimmers,
                                    direction=direction,
                                    max_diff=max_diff,
                                    max_dist=max_dist,
                                    max_repeat=max_repeat)
    return shimmer_alns


def get_shimmer_match_offset(shimmer0, shimmer1,
                             direction=0,
                             max_diff=100,
                             max_dist=1200,
                             max_repeat=1):
    alns = get_shimmer_alns(shimmer0, shimmer1,
                            direction=direction,
                            max_diff=max_diff,
                            max_dist=max_dist,
                            max_repeat=max_repeat)
    if len(alns) == 0:
        return (None, [[]])
    alns.sort(key=lambda x: -len(x[0]))
    aln = alns[0]
    read_offset = aln[0][0][0][3] - aln[0][0][1][3]
    return read_offset, alns


def get_shimmer_match_offset_from_seq(read_seq, ref_seq,
                                      parameters={}):
    alns = get_shimmer_alns_from_seqs(ref_seq, read_seq,
                                      parameters=parameters)
    if len(alns) == 0:
        return (None, [[]])

    alns.sort(key=lambda x: -len(x[0]))
    aln = alns[0]
    read_offset = aln[0][0][0][3] - aln[0][0][1][3]
    return read_offset, alns


def get_align_range(read_seq, ref_seq,
                    read_offset, max_dist=150,
                    aln_str=False):
    read_len = len(read_seq)
    ref_len = len(ref_seq)
    aligned = False
    t_offset = 0
    if aln_str:
        aln_str = 1
    else:
        aln_str = 0
    if read_offset < 0:
        aln = falcon4py.align(read_seq[abs(read_offset):read_len],
                              read_len - abs(read_offset),
                              ref_seq,
                              len(ref_seq),
                              max_dist, aln_str)
        if abs(aln.aln_q_e-aln.aln_q_s) > 500:
            aligned = True
            t_offset = aln.aln_t_s

            # print("t_offset:", t_offset, file=sys.stderr, flush=True)
    else:
        aln = falcon4py.align(read_seq,
                              read_len,
                              ref_seq[read_offset:ref_len],
                              ref_len-read_offset,
                              max_dist, aln_str)
        if abs(aln.aln_q_e-aln.aln_q_s) > 500:
            aligned = True
            t_offset = read_offset
    if not aligned:
        falcon4py.free_alignment(aln)
        aln = None
    return aligned, aln, t_offset


def get_tag_from_seqs(read_seq, ref_seq,
                      read_offset,
                      max_dist=150,
                      aln_len_max_diff=96):
    if read_offset is None:
        return None

    ## print("----",  file=sys.stderr, flush=True)

    aligned, aln, t_offset = get_align_range(read_seq, ref_seq,
                                             read_offset,
                                             max_dist=max_dist,
                                             aln_str=True)
    tag = None
    if aligned:
        dovetail = False
        read_len = len(read_seq)
        ref_len = len(ref_seq)
        diff1 = abs(abs(aln.aln_q_e - aln.aln_q_s) - read_len)
        diff2 = abs(ref_len - read_offset - abs(aln.aln_t_e - aln.aln_t_s))
        if read_offset < 0:
            if abs(abs(aln.aln_q_e - aln.aln_q_s) -
                    (read_len - abs(read_offset))) < aln_len_max_diff:
                dovetail = True
        else:
            if diff1 < aln_len_max_diff or \
               diff2 < aln_len_max_diff:
                dovetail = True

        # print("align range:", read_len, aln.aln_q_s,
        # read_len, aln.aln_q_e, ref_len,
        # aln.aln_t_s, aln.aln_t_e, t_offset, file=sys.stderr, flush=True)
        # print("align range2:", diff1, diff2, file=sys.stderr, flush=True)
        # print("aligned:", read_offset, aligned, dovetail, file=sys.stderr, flush=True)

        if dovetail:
            rng = falcon_ffi.new("aln_range[1]")
            rng[0].s1 = aln.aln_q_s
            rng[0].e1 = aln.aln_q_e
            rng[0].s2 = aln.aln_t_s
            rng[0].e2 = aln.aln_t_e
            tag = falcon4py.get_align_tags(aln.q_aln_str,
                                           aln.t_aln_str,
                                           aln.aln_str_size,
                                           rng, 0, t_offset)

            # print("rng:", aln.aln_q_s, aln.aln_q_e, aln.aln_t_s, aln.aln_t_e, file=sys.stderr, flush=True)

            falcon_ffi.release(rng)
        falcon4py.free_alignment(aln)
    else:
        print("alignment fail in util.get_tag_from_seqs():", read_offset,  file=sys.stderr, flush=True)

    return tag


def get_best_seqs(seqs, best_n=20, levels=2, k=16, w=80):
    all_mmers = []
    mer_count = Counter()
    for i, seq in enumerate(seqs):
        shimmers0 = get_shimmers_from_seq(seq,
                                          rid=0,
                                          levels=levels,
                                          k=k, w=w)
        _shimmers0 = shimmers0.mmers
        seq_mmer_set = set()
        for i in range(_shimmers0.n):
            mmer = mmer2tuple(_shimmers0.a[i])[0]
            seq_mmer_set.add(mmer)
        seq_mmer_set = list(seq_mmer_set)
        mer_count.update(seq_mmer_set)
        all_mmers.append(seq_mmer_set)
        del shimmers0

    match_mmer_count = []
    for i, mmers in enumerate(all_mmers):
        count = 0
        for m in mmers:
            if mer_count[m] >= 2:
                count += 1
        match_mmer_count.append((count, i))

    match_mmer_count.sort(key=lambda x: -x[0])
    seqs_out = []
    for i in range(best_n):
        if i >= len(seqs):
            break
        seqs_out.append(seqs[match_mmer_count[i][1]])

    return seqs_out


def get_cns_from_reads(seqs, sort_reads=True, best_n=20,
                       levels=2, k=16, w=80, max_dist=150,
                       min_cov=1):

    aln_count = 0
    tags = falcon_ffi.new("align_tags_t * [{}]".format(len(seqs)+1))
    if sort_reads:
        seqs = get_best_seqs(seqs, best_n=best_n, levels=levels, k=k, w=w)
    seq0 = seqs[0]
    shimmers0 = get_shimmers_from_seq(seq0,
                                      rid=0,
                                      levels=levels,
                                      k=k, w=w)
    alns = get_shimmer_alns(shimmers0, shimmers0, 0)
    aln = alns[0]
    read_offset = aln[0][0][0][3] - aln[0][0][1][3]
    seq = seq0
    tag = get_tag_from_seqs(seq, seq0, read_offset)
    tags[aln_count] = tag
    aln_count += 1

    for i, seq in enumerate(seqs):
        if i == 0:
            continue

        rid = i * 2
        shimmers_ = get_shimmers_from_seq(seq,
                                          rid=rid,
                                          levels=levels,
                                          k=k, w=w)
        read_offset_n, aln_n = get_shimmer_match_offset(shimmers0, shimmers_)
        del shimmers_

        rseq = rc(seq)
        rid = i * 2 + 1
        shimmers_ = get_shimmers_from_seq(rseq,
                                          rid=rid,
                                          levels=levels,
                                          k=k, w=w)
        read_offset_r, aln_r = get_shimmer_match_offset(shimmers0, shimmers_)
        del shimmers_

        read_offset = 0
        if len(aln_n[0]) == 0 and len(aln_r[0]) == 0:
            continue
        if len(aln_n[0]) > len(aln_r[0]):
            read_offset = read_offset_n
        else:
            read_offset = read_offset_r
            seq = rseq
        if read_offset is None:
            continue
        tag = get_tag_from_seqs(seq, seq0, read_offset,
                                max_dist=max_dist)
        if tag is not None:
            tags[aln_count] = tag
            aln_count += 1

    cns = falcon4py.get_cns_from_align_tags(tags,
                                            aln_count,
                                            len(seq0),
                                            min_cov)
    cns_seq = falcon_ffi.string(cns.sequence)
    falcon4py.free_consensus_data(cns)
    del shimmers0
    for i in range(aln_count):
        falcon4py.free_align_tags(tags[i])
    falcon_ffi.release(tags)

    return cns_seq


SeqIndexData = namedtuple("SeqIndexData", "rname length offset")


class SequenceDatabase(object):
    def __init__(self, index_path="", seqdb_path=""):
        self.basemap = {0: b"N", 1: b"A", 2: b"C", 4: b"G", 8: b"T"}
        self.index_path = index_path
        self.seqdb_path = seqdb_path
        self.name2rid = {}
        self.index_data = {}
        self._load_index()
        self._f = open(seqdb_path, "rb")
        self.seqdb = mmap.mmap(self._f.fileno(), 0,
                               flags=mmap.MAP_SHARED,
                               prot=mmap.PROT_READ)

    def _load_index(self):
        with open(self.index_path) as f:
            for row in f:
                row = row.strip().split()
                rid, rname, rlen, offset = row
                rid = int(rid)
                rlen = int(rlen)
                offset = int(offset)
                self.index_data.setdefault(rid, {})
                self.name2rid[rname] = rid
                self.index_data[rid] = SeqIndexData(rname=rname,
                                                    length=rlen,
                                                    offset=offset)

    def get_subseq_by_rid(self, rid, start=-1, end=-1, direction=0):
        if start == -1 and end == -1:
            start = 0
            end = self.index_data[rid].length
        offset = self.index_data[rid].offset
        s = offset + start
        e = offset + end
        if direction == 1:
            seq = b"".join([self.basemap[(c & 0xF0) >> 4]
                            for c in self.seqdb[s:e]])
        else:
            seq = b"".join([self.basemap[c & 0x0F] for c in self.seqdb[s:e]])
        return seq

    def get_subseq_by_name(self, rname, start=-1, end=-1, direction=0):
        rid = self.name2rid[rname]
        return self.get_subseq_by_rid(rid, start, end, direction=direction)

    def get_seq_index_by_name(self, rname):
        return self.index_data[self.name2rid[rname]]

    def __del__(self):
        self.seqdb.close()
        self._f.close()


