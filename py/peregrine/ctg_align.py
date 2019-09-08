from .utils import get_shimmers_from_seq
from .utils import get_shimmer_alns
from .utils import get_cigar
import vcfpy
from intervaltree import Interval, IntervalTree
from collections import Counter


def get_shimmer_alns_from_seqs(seq0, seq1, parameters={}):
    reduction_factor = parameters.get("reduction_factor", 12)
    direction = parameters.get("direction", 0)
    max_diff = parameters.get("max_diff", 1000)
    max_dist = parameters.get("max_dist", 15000)
    max_repeat = parameters.get("max_repeat", 1)

    seq0_shimmers = get_shimmers_from_seq(
        seq0,
        rid=0,
        reduction_factor=reduction_factor)

    seq1_shimmers = get_shimmers_from_seq(
        seq1,
        rid=1,
        reduction_factor=reduction_factor)

    shimmer_alns = get_shimmer_alns(seq0_shimmers,
                                    seq1_shimmers,
                                    direction=direction,
                                    max_diff=max_diff,
                                    max_dist=max_dist,
                                    max_repeat=max_repeat)
    return shimmer_alns


class SeqDBAligner(object):

    def __init__(self, seqdb0, seqdb1):
        self.sdb0 = seqdb0
        self.sdb1 = seqdb1
        self.default_shimmer_align_parameters = \
            {"reduction_factor": 12,
             "drection": 0,
             "max_diff": 1000,
             "max_dist": 15000,
             "max_repeat": 1}
        self.map_itrees = None

    def _get_vcf_from_cigar(self, seq0, seq1, sname0, sname1, bgn0, bgn1,
                            rpos, qpos, cigars):
        vcf_records = []
        for c in cigars:
            if c[0] == "M":
                rp = rpos
                qp = qpos
                count = 0
                for b0, b1 in zip(seq0[rpos:rpos+c[1]],
                                  seq1[qpos:qpos+c[1]]):
                    if b0 != b1:
                        r = vcfpy.Record(CHROM=sname0, POS=rp + bgn0 + 1,
                                         ID=["."], REF=chr(b0),
                                         ALT=[vcfpy.Substitution("SNV",
                                                                 chr(b1))],
                                         QUAL=50, FILTER=["PASS"],
                                         INFO={"QPOS": qp + bgn1 + 1,
                                               "QNAME": sname1},
                                         FORMAT=["GT"],
                                         calls=[vcfpy.Call("sample0",{"GT":"1|1"})])
                        vcf_records.append(r)
                        count += 1
                    rp += 1
                    qp += 1
                rpos += c[1]
                qpos += c[1]
            elif c[0] == "I":
                if qpos != 0:
                    ins_seq = seq1[qpos-1:qpos+c[1]].decode("ascii")
                else:
                    ins_seq = seq1[qpos:qpos+c[1]].decode("ascii")
                if rpos == 0:
                    r_base = "."
                else:
                    r_base = chr(seq0[rpos-1])

                r = vcfpy.Record(CHROM=sname0, POS=rpos + bgn0, ID=["."],
                                 REF=r_base,
                                 ALT=[vcfpy.Substitution("INS", ins_seq)],
                                 QUAL=50, FILTER=["PASS"],
                                 INFO={"QPOS": qpos + bgn1 + 1,
                                       "QNAME": sname1},
                                        FORMAT=["GT"],
                                        calls=[vcfpy.Call("sample0",{"GT":"1|1"})])
                vcf_records.append(r)
                qpos += c[1]
            elif c[0] == "D":
                if rpos != 0:
                    del_seq = seq0[rpos-1:rpos+c[1]].decode("ascii")
                else:
                    del_seq = seq0[rpos:rpos+c[1]].decode("ascii")
                if qpos == 0:
                    q_base = "."
                else:
                    q_base = chr(seq1[qpos-1])
                r = vcfpy.Record(CHROM=sname0, POS=rpos + bgn0, ID=["."],
                                 REF=del_seq,
                                 ALT=[vcfpy.Substitution("DEL", q_base)],
                                 QUAL=50, FILTER=["PASS"],
                                 INFO={"QPOS": qpos + bgn1 + 1,
                                       "QNAME": sname1},
                                       FORMAT=["GT"],
                                       calls=[vcfpy.Call("sample0",{"GT":"1|1"})])
                vcf_records.append(r)
                rpos += c[1]
        return vcf_records

    def get_vcf_records(self, seq0_info, seq1_info,
                        parameters=None, extend=True, ext_size=5000):
        if parameters is None:
            parameters = self.default_shimmer_align_parameters
        direction = parameters.get("direction", 0)
        score = parameters.get("score", (2, -4, 4, 2))
        ctg_direction = direction
        sname0, bgn0, end0 = seq0_info
        sname1, bgn1, end1 = seq1_info
        seq0 = self.sdb0.get_subseq_by_name(sname0, bgn0, end0)
        seq1 = self.sdb1.get_subseq_by_name(sname1, bgn1, end1,
                                            direction=direction)
        if bgn0 == -1:
            bgn0 = 0
        if bgn1 == -1:
            bgn1 = 0
        if end0 == -1:
            end0 = self.sdb0.get_seq_index_by_name(sname0).length
        if end1 == -1:
            end1 = self.sdb0.get_seq_index_by_name(sname1).length

        if direction == 1:
            parameters["direction"] = 0

        shimmer_alns = get_shimmer_alns_from_seqs(seq0, seq1, parameters)

        vcf_records = []
        for aln, aln_d in shimmer_alns:
            if len(aln) < 2:
                continue
            if extend is True and len(vcf_records) == 0:
                ext = min(ext_size, aln[0][0].pos_end,  aln[0][1].pos_end)
                x0, x1 = aln[0][0].pos_end - ext, aln[0][0].pos_end
                if x0 < 0:
                    x0 = 0
                y0, y1 = aln[0][1].pos_end - ext, aln[0][1].pos_end
                if y0 < 0:
                    y0 = 0
                rpos = x0
                qpos = y0
                sub_seq0 = seq0[x0:x1]
                sub_seq1 = seq1[y0:y1]
                cigars, aln_score = get_cigar(sub_seq0, sub_seq1,
                                              score=score)
                v = self._get_vcf_from_cigar(seq0, seq1, sname0, sname1,
                                             bgn0, bgn1,
                                             rpos, qpos, cigars)
                vcf_records.append((((x0, x1), (y0, y1), 0), v, cigars))

            for i in range(len(aln)-1):
                x0, x1 = aln[i][0].pos_end, aln[i+1][0].pos_end
                y0, y1 = aln[i][1].pos_end, aln[i+1][1].pos_end
                aln_range = ((x0, x1),
                             (y0, y1), len(aln))
                sub_seq0 = seq0[x0:x1]
                sub_seq1 = seq1[y0:y1]
                cigars, aln_score = get_cigar(sub_seq0, sub_seq1,
                                              score=score)
                v = self._get_vcf_from_cigar(seq0, seq1, sname0, sname1,
                                             bgn0, bgn1,
                                             x0, y0, cigars)
                vcf_records.append((aln_range, v, cigars))

        if len(vcf_records) > 0 and extend is True:
            ext = min(ext_size, end0-x1, end1-y1)
            x0, x1 = x1, x1 + ext
            y0, y1 = y1, y0 + y1 + ext
            slen = self.sdb1.get_seq_index_by_name(sname1).length
            if y1 > slen:
                y1 = slen
            rpos = x0
            qpos = y0
            sub_seq0 = seq0[x0:x1]
            sub_seq1 = seq1[y0:y1]
            cigars, aln_score = get_cigar(sub_seq0, sub_seq1,
                                          score=score)
            v = self._get_vcf_from_cigar(seq0, seq1, sname0, sname1,
                                         bgn0, bgn1,
                                         rpos, qpos, cigars)
            vcf_records.append((((x0, x1), (y0, y1), 0), v, cigars))

        return vcf_records, ctg_direction

    def load_map_file(self, map_file_path):
        self.map_itrees = {}
        with open(map_file_path) as f:
            for row in f:
                row = row.strip().split()
                (ref_id, ref_bgn, ref_end,
                    ctg_id, ctg_bgn, ctg_end,
                    ctg_direction, mcount0, mcount1) = [int(_) for _ in row]
                self.map_itrees.setdefault(ref_id, IntervalTree())
                self.map_itrees[ref_id][ref_bgn:ref_end] = \
                    (ctg_id, ctg_bgn, ctg_end,
                     ctg_direction, mcount0, mcount1)

    def _generate_interval_candidate(self, shimmer_alns,
                                     seq0, seq1,
                                     itvl,
                                     bgn0_, bgn0, end0,
                                     min_size=500):

        candidate = [None, None]
        for aln, aln_dist in shimmer_alns:

            x0 = aln[0].mmer0.pos_end
            # y0 = aln[0].mmer1.pos_end
            x1 = aln[-1].mmer0.pos_end
            # y1 = aln[-1].mmer1.pos_end

            if (x1-x0) < min_size:
                continue

            for i in range(len(aln)-1):
                x0 = aln[i].mmer0.pos_end
                y0 = aln[i].mmer1.pos_end
                x1 = aln[i+1].mmer0.pos_end
                y1 = aln[i+1].mmer1.pos_end

                if (bgn0 < x0 + bgn0_ or x1 + bgn0_ < bgn0) and \
                   (end0 < x0 + bgn0_ or x1 + bgn0_ < end0):
                    continue

                cigars, aln_score = get_cigar(seq0[x0:x1], seq1[y0:y1])

                rpos = x0 + bgn0_
                qpos = y0 + itvl.begin

                for cigar in cigars:
                    if cigar[0] in ('M', 'I'):
                        if rpos < bgn0 and bgn0 < rpos + cigar[1]:
                            pos = qpos + bgn0 - rpos
                            candidate[0] = pos
                        if rpos < end0 and end0 < rpos + cigar[1]:
                            pos = qpos + end0 - rpos
                            candidate[1] = pos
                        rpos += cigar[1]
                    if cigar[0] in ('M', 'D'):
                        qpos += cigar[1]
        return candidate

    def get_target_itree(self, intervals, sid=0, padding=5000):
        target_itrees = {}
        for interval in sorted(intervals):
            (ctg_id, ctg_bgn, ctg_end,
                ctg_direction, mcount0, mcount1) = interval.data
            target_itrees.setdefault(ctg_id, IntervalTree())
            b = ctg_bgn - padding
            e = ctg_bgn + padding
            target_itrees[ctg_id][b:e] = \
                (sid, ctg_direction, interval.begin, interval.end)
        for t_id in target_itrees:
            target_itrees[t_id].merge_overlaps(
                data_reducer=lambda x, y: x+[y[1]],
                data_initializer=[],
                strict=False)
        return target_itrees

    def map_small_interval(self, seq0_info, padding=5000):
        """
        find interval in seq1 that is corresponding to
        sname0:bgn0-bgn1
        bgn1 - bgn0 should be less than 250000 for now
        """
        assert self.map_itrees is not None
        sname0, bgn0, end0 = seq0_info
        assert end0 - bgn0 > 0
        assert end0 - bgn0 < 250000
        sid = self.sdb0.name2rid[sname0]

        bgn0_ = bgn0 - padding
        end0_ = end0 + padding
        if bgn0_ < 0:
            bgn0_ = 0
        if end0_ > self.sdb0.index_data[sid].length:
            end0_ = self.sdb0.index_data[sid].length
        seq0 = self.sdb0.get_subseq_by_rid(sid, bgn0_, end0_)

        target_itrees = self.get_target_itree(
                                self.map_itrees[sid][bgn0_:end0_],
                                sid, padding)

        candidates = []
        for t_id in target_itrees:
            for itvl in sorted(target_itrees[t_id]):
                direction = Counter(itvl.data).most_common(1)[0][0]
                seq1 = self.sdb1.get_subseq_by_rid(t_id,
                                                   itvl.begin, itvl.end,
                                                   direction=direction)
                shimmer_alns = get_shimmer_alns_from_seqs(
                                    seq0, seq1,
                                    parameters={"reduction_factor": 4})
                pos1, pos2 = self._generate_interval_candidate(
                                      shimmer_alns,
                                      seq0, seq1,
                                      itvl,
                                      bgn0_, bgn0, end0,
                                      min_size=2*padding)
                if pos1 is not None and pos2 is not None:
                    if direction == 1:
                        pos1 = self.sdb1.index_data[t_id].length - pos1
                        pos2 = self.sdb1.index_data[t_id].length - pos2

                    candidate = [self.sdb1.index_data[t_id].rname,
                                 direction, pos1, pos2]
                    candidates.append(candidate)

        return candidates
