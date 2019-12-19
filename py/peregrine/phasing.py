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
        for i in range(len(self.seqs)):
            rid = i
            seq = self.seqs[i]
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

    def get_read_group_graph(self):
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
                if score > 0.99 and len(d_ovlp) > 3:
                    G.add_edge(r0, r1, score=score)
        return G, scores
