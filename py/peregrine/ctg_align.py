from .utils import get_shimmers_from_seq
from .utils import get_shimmer_alns
from .utils import get_cigar
import vcfpy


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

    def get_shimmer_alns(self, seq0, seq1, parameters=None):
        if parameters is None:
            parameters = self.default_shimmer_align_parameters
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

    def _get_vcf_from_cigar(self, seq0, seq1, sname0, bgn0, bgn1,
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
                                            INFO={"QPOS": [qp + bgn1 + 2]})
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
                                 INFO={"QPOS": [qpos + bgn1 + 2]})
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
                                 INFO={"QPOS": [qpos + bgn1 + 2]})
                vcf_records.append(r)
                rpos += c[1]
        return vcf_records

    def get_vcf_records(self, seq0_info, seq1_info,
                        parameters=None, extend=True):
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
        if direction == 1:
            parameters["direction"] = 0

        shimmer_alns = self.get_shimmer_alns(seq0, seq1, parameters)

        vcf_records = []
        for aln, aln_d in shimmer_alns:
            if len(aln) < 2:
                continue
            if extend is True and len(vcf_records) == 0:
                x0, x1 = 0, aln[0][0].pos_end
                y1 = aln[0][1].pos_end
                y0 = y1 - x1 - 100
                if y0 < 0:
                    y0 = 0
                rpos = x0
                qpos = y0
                sub_seq0 = seq0[x0:x1]
                sub_seq1 = seq1[y0:y1]
                cigars, aln_score = get_cigar(sub_seq0, sub_seq1,
                                              score=score)
                v = self._get_vcf_from_cigar(seq0, seq1, sname0,
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
                v = self._get_vcf_from_cigar(seq0, seq1, sname0,
                                                bgn0, bgn1,
                                                x0, y0, cigars)
                vcf_records.append((aln_range, v, cigars))

        if extend is True:
            x0, x1 = x1, end0 - bgn0
            y0 = y1
            y1 = y0 + x1 - x0 + 100
            slen = self.sdb1.get_seq_index_by_name(sname1).length
            if y1 > slen:
                y1 = slen
            rpos = x0
            qpos = y0
            sub_seq0 = seq0[x0:x1]
            sub_seq1 = seq1[y0:y1]
            cigars, aln_score = get_cigar(sub_seq0, sub_seq1,
                               score=score)
            v = self._get_vcf_from_cigar(seq0, seq1, sname0,
                                            bgn0, bgn1,
                                            rpos, qpos, cigars)
            vcf_records.append((((x0, x1), (y0, y1), 0), v, cigars))

        return vcf_records, ctg_direction
