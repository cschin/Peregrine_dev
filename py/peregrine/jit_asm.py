from .utils import get_shimmers_from_seq
from .utils import rc
from .utils import mmer2tuple
# from .utils import SequenceDatabase
# from .utils import get_shimmer_alns
# from .utils import get_cigar
# from .utils import get_shimmer_alns_from_seqs

import networkx as nx
from collections import Counter
from networkx.algorithms.dag import transitive_reduction
from networkx.algorithms.dag import dag_longest_path
from networkx import weakly_connected_components


class JITAssembler(object):
    def __init__(self, reads, k=16, w=80, min_shimmer_dist=250):
        """
        reads should be a directionary mapping
        interger reads to byte read seqeunces
        """
        self.reads = reads
        self.k = k
        self.w = w
        self.min_shimmer_dist = min_shimmer_dist
        self.read2smps = {}
        self.smp2reads = {}
        self.read2shimmers = {}
        self.smp_count = None
        self.full_graph = None
        self.circles = None
        self.DAG = None
        self.anchor2read = None

    def _get_read_shimmer_maps(self,
                               min_dist=250,
                               remove_repeat=False):
        """
        Loop through the reads, create the shimmer for each reads
        and build the map between the shimmers to reads and the
        reads to shimmers
        """
        for rid, seq in self.reads.items():
            shimmers0 = get_shimmers_from_seq(seq,
                                              rid=rid,
                                              levels=2,
                                              k=self.k,
                                              w=self.w)
            self.read2shimmers[rid] = shimmers0
            mmpairs0 = []
            mmpairs0_pos = {}
            mmpairs1 = []
            mmpairs1_pos = {}
            for i in range(shimmers0.mmers.n-1):
                mmer0 = mmer2tuple(shimmers0.mmers.a[i])
                mmer1 = mmer2tuple(shimmers0.mmers.a[i+1])
                if abs(mmer0.pos_end - mmer1.pos_end) < min_dist:
                    continue
                if mmer0.mmer == mmer1.mmer:
                    continue
                key = (mmer0.mmer, mmer1.mmer)
                pos = (mmer0.pos_end, mmer1.pos_end)
                mmpairs0.append(key)
                mmpairs0_pos.setdefault(key, [])
                mmpairs0_pos[key].append(pos)

                key = (mmer1.mmer, mmer0.mmer)
                pos = (mmer1.pos_end - self.k, mmer0.pos_end - self.k)
                mmpairs1.append(key)
                mmpairs1_pos.setdefault(key, [])
                mmpairs1_pos[key].append(pos)

            self.read2smps[(rid, 0)] = []
            self.read2smps[(rid, 1)] = []

            for key in mmpairs0:
                pos = mmpairs0_pos[key]
                # may filter out duplicated SMPs in a read
                if remove_repeat and len(pos) > 1:
                    continue
                pos = pos[0]
                self.smp2reads.setdefault(key, [])
                self.smp2reads[key].append((rid, 0, pos))
                self.read2smps[(rid, 0)].append(key)

            for key in mmpairs1[::-1]:
                pos = mmpairs1_pos[key]
                # may filter out duplicated SMPs in a read
                if remove_repeat and len(pos) > 1:
                    continue
                pos = pos[0]
                self.smp2reads.setdefault(key, [])
                self.smp2reads[key].append((rid, 1, pos))
                self.read2smps[(rid, 1)].append(key)

        self.smp_count = Counter()
        for v in self.read2smps.values():
            self.smp_count.update(v)

    def _get_quick_anchor_graph(self,
                                min_smp_count=8,
                                max_node_count=80,
                                min_edge_count=4):
        """
        Create the interlocked draft anchor graph
        An "anchor" is a shimmer pair which is near
        the end of a read but also as an internal
        shimmer pair of other reads
        """
        for rid, smps in self.read2smps.items():
            self.read2smps[rid] = [s for s in smps
                                   if self.smp_count[s] >= min_smp_count]

        inner_smps = set()

        for rid, smps in self.read2smps.items():
            if len(smps) < 3:
                continue
            for i in range(1, len(smps)-1):
                inner_smps.add(smps[i])

        anchoring_smps = {}
        self.anchor2read = {}
        for rid, smps in self.read2smps.items():
            if len(smps) < 2:
                continue
            if smps[0] in inner_smps and smps[-1] in inner_smps:
                anchoring_smps.setdefault(smps[0], [])
                anchoring_smps[smps[0]].append(rid)
                anchoring_smps.setdefault(smps[-1], [])
                anchoring_smps[smps[-1]].append(rid)
                self.anchor2read.setdefault(smps[0], [])
                self.anchor2read[smps[0]].append((rid, "B"))
                self.anchor2read.setdefault(smps[-1], [])
                self.anchor2read[smps[-1]].append((rid, "E"))
        node_count = Counter()
        edges = []
        for rid, smps in self.read2smps.items():
            path = []
            for smp in smps:
                if smp in anchoring_smps:
                    path.append(smp)
            node_count.update(path)

        for rid, smps in self.read2smps.items():
            path = []
            for smp in smps:
                if smp in anchoring_smps and node_count[smp] <= max_node_count:
                    path.append(smp)
            for i in range(len(path)-1):
                edges.append((path[i], path[i+1]))

        edge_count = Counter(edges)
        G = nx.DiGraph()
        for rid, smp in self.read2smps.items():
            path = []
            for smp in smp:
                if smp in anchoring_smps:
                    path.append(smp)
            for i in range(len(path)-1):
                edge = (path[i], path[i+1])
                if edge_count[edge] < min_edge_count:
                    continue
                if G.has_edge(*edge):
                    G.edges[edge]["count"] += 1
                else:
                    G.add_edge(*edge, count=1, in_cycle=0)

        self.circles = []
        for e in G.edges():
            G.edges[e]["in_cycle"] = 0
        G_ = G.copy()
        while 1:
            try:
                C = nx.find_cycle(G_)
                self.circles.append(C)
                for edge in C:
                    G.edges[edge]["in_cycle"] = 1
                    G_.remove_edge(*edge)
            except nx.NetworkXNoCycle:
                break

        self.full_graph = G.copy()
        self.DAG = transitive_reduction(G_)

    def _get_ctgs_from_DAG(self):
        def reverse_smp(smp):
            rsmp = tuple(list(smp)[::-1])
            return rsmp

        span_read_count = Counter()

        for subG in [self.DAG.subgraph(c) for c in weakly_connected_components(self.DAG)]:
            path = dag_longest_path(subG)
            v = path[0]
            for w in path[1:]:
                v_reads = set([(x[0], x[1]) for x in self.smp2reads[v]])
                w_reads = set([(x[0], x[1]) for x in self.smp2reads[w]])
                span_read_count.update(v_reads & w_reads)
                v = w

        used_edges = set()
        ctgs = []
        edge_count = 0
        for subG in [self.DAG.subgraph(c) for c in weakly_connected_components(self.DAG)]:
            subG = subG.copy()
            while len(subG.edges) != 0:
                path = dag_longest_path(subG)
                v = path[0]
                seqs = []
                added_edges = []
                for w in path[1:]:
                    edge_count += 1
                    if (v, w) in used_edges:
                        added_edges.append((v, w))
                        if len(seqs) > 1:
                            ctgs.append(b"".join(seqs))
                        seqs = []
                        continue

                    v_reads = set([x[0:2] for x in self.smp2reads[v]])
                    w_reads = set([x[0:2] for x in self.smp2reads[w]])

                    count_reads = []
                    for n in list(v_reads & w_reads):
                        count_reads.append((span_read_count[n], n))
                    count_reads.sort()

                    if len(count_reads) == 0:
                        if len(seqs) > 1:
                            ctgs.append(b"".join(seqs))
                        seqs = []
                        continue

                    m0 = [x for x in self.smp2reads[v]
                          if x[0:2] == count_reads[-1][1]][0]
                    m1 = [x for x in self.smp2reads[w]
                          if x[0:2] == count_reads[-1][1]][0]
                    read_id = m0[0]
                    direction = m0[1]
                    pos0 = m0[2][0]
                    pos1 = m1[2][0]
                    if direction == 1:
                        seq = rc(self.reads[read_id][pos1:pos0])
                    else:
                        seq = self.reads[read_id][pos0:pos1]
                    seqs.append(seq)

                    used_edges.add((v, w))
                    used_edges.add((reverse_smp(w), reverse_smp(v)))
                    added_edges.append((v, w))
                    v = w

                if len(seqs) > 0:
                    ctgs.append(b"".join(seqs))

                for v, w in added_edges:
                    if subG.has_edge(v, w):
                        subG.remove_edge(v, w)
        return ctgs
