from typing import List
from .utils import rc
from .utils import get_shimmer_match_offset_from_seq
from .utils import get_tag_from_seqs
from collections import Counter
import networkx as nx


class GraphPhaser(object):
    """
    A simple Graph bases phasing algorithm
    supporting multiple haplotype (>=2) phasing.
    It is mostly experimental and based on
    huristics for now.
    """
    def __init__(self, seqs: List[bytes],
                 ref_seq: bytes,
                 min_count=4,
                 max_off=-1):
        """
        seqs: a list of byte sequence
        """
        self.seqs = seqs
        self.ref_seq = ref_seq
        self.delta = []
        self.tpos2rid = {}
        self.tags = []
        self.delta2read = {}
        self.read2delta = {}
        self.delta2reads2 = {}
        self.read2deltas2 = {}
        self.min_count = min_count
        self.max_off = max_off

    def _filter_deltas(self):
        min_count = self.min_count
        max_off = self.max_off
        self.deltas.sort(key=lambda x: x[1][0])
        delta2reads = {}
        for rid, d in self.deltas:
            delta2reads.setdefault(d, set())
            delta2reads[d].add(rid)

        delta2reads2 = {}
        read2deltas2 = {}

        for d, reads in sorted(delta2reads.items()):
            t_pos = d[0]
            if t_pos not in self.tpos2rid:
                continue

            c0, c1 = len(self.tpos2rid[t_pos]), len(reads)

            if c1 <= min_count or c0 - c1 <= min_count:
                continue

            if max_off > 0:
                if abs(c1/c0 - 0.5) > max_off:
                    continue

            if d[1] <= 1:
                delta2reads2[(d, 1)] = set(reads)
                delta2reads2[(d, 0)] = set(self.tpos2rid[t_pos])-set(reads)
                for r in list(reads):
                    read2deltas2.setdefault(r, set())
                    read2deltas2[r].add((d, 1))
                for r in list(set(self.tpos2rid[t_pos])-set(reads)):
                    read2deltas2.setdefault(r, set())
                    read2deltas2[r].add((d, 0))
        self.delta2reads = delta2reads
        self.delta2reads2 = delta2reads2
        self.read2deltas2 = read2deltas2

    def _get_deltas(self):

        deltas = []
        tpos2rid = {}
        tags = []

        ref_seq = self.ref_seq

        reads = set()
        for rid in self.seqs:
            seq = self.seqs[rid]
            rseq = rc(seq)
            read_offset0, aln0 = \
                get_shimmer_match_offset_from_seq(seq, ref_seq,
                                                  parameters={
                                                  "reduction_factor": 3})
            read_offset1, aln1 = \
                get_shimmer_match_offset_from_seq(rseq, ref_seq,
                                                  parameters={
                                                  "reduction_factor": 3})
            if len(aln0[0]) > len(aln1[0]):
                read_offset = read_offset0
            else:
                read_offset = read_offset1
                seq = rseq
            if read_offset is None:
                continue
            tag = get_tag_from_seqs(seq,
                                    ref_seq,
                                    read_offset,
                                    max_dist=150,
                                    aln_len_max_diff=120)
            if tag is None:
                continue
            for i in range(tag.len):
                ctag = tag.align_tags[i]
                if ctag.delta == 0:
                    tpos2rid.setdefault(ctag.t_pos, set())
                    tpos2rid[ctag.t_pos].add(rid)
                if ctag.delta > 0 or ctag.q_base == b"-":
                    deltas.append((rid, (ctag.t_pos, ctag.delta, ctag.q_base)))
                    reads.add(rid)
        self.deltas = deltas
        self.tpos2rid = tpos2rid
        self.tags = tags
        self._filter_deltas()
        return deltas, reads

    def get_read_group_graph(self, min_score=0.99, min_overlap=3):
        self._get_deltas()

        read2start = {}
        read2end = {}
        read_deltas = []
        for r, d in self.read2deltas2.items():
            d = sorted(list(d))
            # print(d)
            read2start[r] = d[0][0]
            read2end[r] = d[-1][0]
            read_deltas.append((r, d))
        read_deltas.sort(key=lambda x: x[1][0])

        G = nx.DiGraph()
        scores = []
        for i in range(len(read_deltas)):
            for j in range(i+1, len(read_deltas)):
                r0, d0 = read_deltas[i]
                r1, d1 = read_deltas[j]
                if read2start[r1] > read2end[r0]:
                    break
                d_ovlp = {}
                for _ in d0:
                    if _[0] < read2start[r1] or _[0] > read2end[r1]:
                        continue
                    d_ovlp.setdefault(_, [0, 0])
                    d_ovlp[_][0] = 1
                for _ in d1:
                    if _[0] < read2start[r0] or _[0] > read2end[r0]:
                        continue
                    d_ovlp.setdefault(_, [0, 0])
                    d_ovlp[_][1] = 1
                count = Counter([tuple(_) for _ in d_ovlp.values()])
                score = count.get((0, 0), 0) + count.get((1, 1), 0)
                score -= count.get((0, 1), 0) + count.get((1, 0), 0)
                scores.append(score/len(d_ovlp))
                score /= len(d_ovlp)
                # print(i, j, score)
                if score >= min_score and len(d_ovlp) >= min_overlap:
                    G.add_edge(r0, r1, score=score)
        return G, scores

    def biallelic_phasing(self, niter=3):

        def get_read_phases(reads, read2phase):
            phases = []
            for r in reads:
                if r in read2phase:
                    phases.append(read2phase[r])
            return phases

        read2phase = {}
        delta2phase = {}
        deltas_ = set([_[0] for _ in self.delta2reads2.keys()])
        deltas_ = list(deltas_)
        deltas_.sort()
        last_bgn_delta = None
        pre_delta = None
        blocks = []
        for i in range(niter):
            for delta in deltas_:
                delta0 = (delta, 0)
                delta1 = (delta, 1)
                reads0 = self.delta2reads2[delta0]
                reads1 = self.delta2reads2[delta1]
                phases0 = get_read_phases(reads0, read2phase)
                phases1 = get_read_phases(reads1, read2phase)

                delta_phase0 = delta2phase.get(delta0, None)
                delta_phase1 = delta2phase.get(delta1, None) 

                # New block
                if delta_phase0 == None and \
                   delta_phase1 == None and \
                   len(phases0) == 0 and \
                   len(phases1) == 0:
                    for r in reads0:
                        read2phase[r] = 0
                    for r in reads1:
                        read2phase[r] = 1
                    phases0 = get_read_phases(reads0, read2phase)
                    phases1 = get_read_phases(reads1, read2phase)
                    if last_bgn_delta is not None:
                        blocks.append((last_bgn_delta[0], pre_delta[0]))
                    last_bgn_delta = delta

                p0count = Counter(phases0)
                p1count = Counter(phases1)

                d0 = p0count.get(0, 0) + p1count.get(1, 0) \
                    - p0count.get(1, 0) - p1count.get(0, 0)

                if d0 > 0:
                    delta2phase[delta0] = 0
                    delta2phase[delta1] = 1
                    for r in reads0:
                        if r not in read2phase:
                            read2phase[r] = 0
                    for r in reads1:
                        if r not in read2phase:
                            read2phase[r] = 1     
                else:
                    delta2phase[delta0] = 1
                    delta2phase[delta1] = 0
                    for r in reads0:
                        if r not in read2phase:
                            read2phase[r] = 1
                    for r in reads1:
                        if r not in read2phase:
                            read2phase[r] = 0
                pre_delta = delta
            if i == 0 and last_bgn_delta is not None and pre_delta is not None:
                blocks.append((last_bgn_delta[0], pre_delta[0]))
            for r in read2phase:
                phases = [delta2phase[_] for _ in list(self.read2deltas2[r])]
                read2phase[r] = Counter(phases).most_common(1)[0][0]

        return blocks, delta2phase, read2phase
