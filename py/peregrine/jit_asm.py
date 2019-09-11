from .utils import SequenceDatabase
from .utils import get_shimmers_from_seq
from .utils import get_shimmer_alns
from .utils import get_cigar
from .utils import get_shimmer_alns_from_seqs
from .utils import rc
from .utils import mmer2tuple
import networkx as nx
from collections import Counter
from networkx.algorithms.dag import transitive_reduction


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

    def _get_read_shimmer_maps(self):
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
            mmpairs1 = []
            for i in range(shimmers0.mmers.n-1):
                mmer0 = mmer2tuple(shimmers0.mmers.a[i])
                mmer1 = mmer2tuple(shimmers0.mmers.a[i+1])
                if abs(mmer0.pos_end - mmer1.pos_end) < 250:
                    continue
                if mmer0.mmer == mmer1.mmer:
                    continue
                key = (mmer0.mmer, mmer1.mmer)
                mmpairs0.append(key)
                self.smp2reads.setdefault(key, [])
                self.smp2reads[key].append((rid, 1))

                key = (mmer1.mmer, mmer0.mmer)
                mmpairs1.append(key)
                self.smp2reads.setdefault(key, [])
                self.smp2reads[key].append(rid*10+1)

            self.read2smps[(rid, 0)] = ["{}:{}".format(*x)
                                        for x in mmpairs0]
            self.read2smps[(rid, 1)] = ["{}:{}".format(*x)
                                        for x in mmpairs1[::-1]]
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
        the end of a read but alsoe as an internal
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
        for v in G.nodes():
            G.nodes[edge[0]]["in_cycle"] = 0
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
