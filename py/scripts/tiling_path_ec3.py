#!/usr/bin/env python3

import mmap
import sys
import numpy
from peregrine._falcon4py import ffi as falcon_ffi
from peregrine._falcon4py import lib as falcon4py
from peregrine._shimmer4py import ffi as shimmer_ffi
from peregrine._shimmer4py import lib as shimmer4py
from peregrine.utils import SequenceDatabase
#from peregrine.utils import rc
from peregrine.utils import get_shimmers_from_seq
from peregrine.utils import get_shimmer_alns
from peregrine.utils import get_tag_from_seqs
#from peregrine.utils import get_shimmer_match_offset
from peregrine.utils import get_shimmer_alns_from_seqs
#from peregrine._shimmer4py import ffi, lib

basemap = {1:"A",2:"C",4:"G",8:"T"}
stitching_overhang_size = 500

def get_cns(template_seq, supporting_reads, max_dist = 128, min_cov=1):

    aln_count = 0
    tags = falcon_ffi.new("align_tags_t * [{}]".format(len(supporting_reads)+1))
    seq0 = template_seq
    tag = get_tag_from_seqs(seq0, seq0, 0)
    tags[aln_count] = tag
    aln_count += 1

    for rid, read_offset, seq in supporting_reads:

        tag = get_tag_from_seqs(seq, seq0, read_offset,
                                max_dist=max_dist)
        if tag is not None:
            tags[aln_count] = tag
            aln_count += 1

    print("start generating cns", file=sys.stderr, flush=True)
    cns = falcon4py.get_cns_from_align_tags(tags,
                                            aln_count,
                                            len(seq0),
                                            min_cov)
    print("finish generating cns", file=sys.stderr, flush=True)

    cns_seq = falcon_ffi.string(cns.sequence)
    falcon4py.free_consensus_data(cns)
    for i in range(aln_count):
        falcon4py.free_align_tags(tags[i])
    falcon_ffi.release(tags)

    return cns_seq

def get_offset(template_mmers, read_seq, l, w, k):
    mmers = shimmer_ffi.new("mm128_v *")

    shimmer4py.mm_sketch(c_null, read_seq, l, w, k, 0, 0, mmers)
    aln = shimmer4py.shmr_aln(template_mmers, mmers, 0, 100, 1200, 1)

    max_chain_idx = -1
    max_chain_n = -1
    #print("X:", aln.n)
    if aln.n == 0:
        shimmer4py.free(mmers.a)
        shimmer_ffi.release(mmers)
        shimmer4py.free_shmr_alns(aln)
        return None
    for idx in range(aln.n):
        if aln.a[idx].idx0.n > max_chain_n:
            max_chain_n = aln.a[idx].idx0.n
            max_chain_idx = idx
    idx0 = aln.a[max_chain_idx].idx0.a[0] 
    idx1 = aln.a[max_chain_idx].idx1.a[0]
    y = template_mmers.a[idx0].y
    template_pos_end = ((y & 0xFFFFFFFF) >> 1) + 1
    y = mmers.a[idx1].y
    read_pos_end = ((y & 0xFFFFFFFF) >> 1) + 1
    shimmer4py.free(mmers.a)
    shimmer_ffi.release(mmers)
    shimmer4py.free_shmr_alns(aln)
    offset = template_pos_end - read_pos_end
    return offset

def stiching_reads(tiling_path_data, seqdb, ovlps):
    segments = []
    # I don't like to have the first read as it breaks the string formulation,
    # but poeple like it for no reason, so I will just do it
    ctg_id, v, w, r, s, e, olen, idt, _1, _2 = tiling_path_data[0]
    v = v.split(":")
    rid0 = int(v[0])
    s0 = seqdb.index_data[rid0].offset 
    slen0 = seqdb.index_data[rid0].length
    e0 = s0 + slen0
    bseq0 = seqdb.seqdb[s0:e0]
    strand0 = 0 if v[1] == "E" else 1

    seq = shimmer_ffi.new("char[{}]".format(slen0))
    shimmer4py.decode_biseq(bseq0,
            seq,
            slen0,
            strand0)

    ctg_len = len(seq)
    segments.append((ctg_len, 0, seq))
    print("start stiching, ctg:", ctg_id, file=sys.stderr, flush=True)

    supporting_seq_ids = set()
    for row in tiling_path_data:
        print(" ".join(row), file=sys.stderr, flush=True)
        ctg_id, v, w, r, s, e, olen, idt, _1, _2 = row
        v = v.split(":")
        w = w.split(":")
        s = int(s)
        e = int(e)
        olen = int(olen)
        idt = float(idt)

        rid0 = int(v[0])

        tiling_read_direction = 0 if int(e) > int(s) else 1
        for rr in ovlps[rid0]:
            _, _, s_, e_, direction = ovlps[rid0][rr]
            if abs(s_-e_) < 4000:
                continue
            #direction = ovlps[rid0][rr][-1]
            supporting_seq_ids.add( (rr, tiling_read_direction, direction) ) 

        s0 = seqdb.index_data[rid0].offset 
        slen0 = seqdb.index_data[rid0].length
        e0 = s0 + slen0
        bseq0 = seqdb.seqdb[s0:e0]
        strand0 = 0 if v[1] == "E" else 1

        rid1 = int(w[0])
        s1 = seqdb.index_data[rid1].offset 
        slen1 = seqdb.index_data[rid1].length
        e1 = s1 + slen1
        bseq1 = seqdb.seqdb[s1:e1]
        strand1 = 0 if w[1] == "E" else 1

        offset1 = slen0 - stitching_overhang_size
        offset2 = slen1 - abs(e-s) - stitching_overhang_size
        match = shimmer4py.ovlp_match(bseq0[offset1:], slen0 - offset1, strand0,
                               bseq1[offset2:], slen1 - offset2, strand1,
                               100)

        if strand1 == 1:
            s, e = slen1 - s, slen1 - e
        assert(e > s)
        seg_size = e - s + stitching_overhang_size - match.t_m_end
        seq = shimmer_ffi.new("char[{}]".format(seg_size))
        shimmer4py.decode_biseq(bseq1[e-seg_size:e],
                seq,
                seg_size,
                strand1)
        segments.append((ctg_len,
            ctg_len - stitching_overhang_size + match.q_m_end,
            seq))
        ctg_len -= (stitching_overhang_size - match.q_m_end)
        ctg_len += (stitching_overhang_size - match.t_m_end) + e - s

        shimmer4py.free_ovlp_match(match)

    ctg_str = numpy.ones(ctg_len, dtype=numpy.byte)
    ctg_str *= ord('N')
    for seg in segments:
        s = seg[1]
        e = seg[1] + len(shimmer_ffi.string(seg[2]))
        ctg_str[s:e] = list(shimmer_ffi.string(seg[2]))
        shimmer_ffi.release(seg[2])
    template_seq = bytes(ctg_str)
    print("stiching done, ctg:", ctg_id, file=sys.stderr, flush=True)
    #print("template seq:", template_seq, file=sys.stderr, flush=True)

    return template_seq, supporting_seq_ids


def get_ovlps(overlap_data_fn, tiling_reads):
    ovlps = {}
    with open(overlap_data_fn) as f: 
        for r in f:
            r = r.strip().split()
            if r[0] == "-":
                continue
            id1, id2 = r[:2]
            id1 = int(id1)
            id2 = int(id2)
            
            if id1 in tiling_reads:
                ref_s = int(r[5])
                ref_e = int(r[6])
                direction = int(r[8])
                s = int(r[9])
                e = int(r[10])
                l = int(r[11])
                if direction == 1:
                    s, e = l - e, l - s
                ovlps.setdefault(id1, {})
                ovlps[id1][id2] = (ref_s, ref_e, s, e, direction)

            if id2 in tiling_reads:
                ref_s = int(r[9])
                ref_e = int(r[10])
                # ref_l = int(r[11])
                direction = int(r[8])
                s = int(r[5])
                e = int(r[6])
                l = int(r[7])
                if direction == 1:
                    s, e = l - e, l - s
                ovlps.setdefault(id2, {})
                ovlps[id2][id1] = (ref_s, ref_e, s, e, direction)
    return ovlps


def get_tiling_reads(tiling_path_fn, total_chunk, my_chunk):
    tiling_path_data = {}
    tiling_reads = set()
    with open(tiling_path_fn) as f:
        for row in f:
            row = row.strip().split()
            ctg = row[0]
            ctg_hash = int(ctg[:6]) 
            if ctg_hash % total_chunk != my_chunk % total_chunk:
                continue
            tiling_path_data.setdefault(row[0], [])
            tiling_path_data[row[0]].append(row)
            ctg_id, v, w, r, s, e, olen, idt, _1, _2 = row
            v = v.split(":")[0]
            w = w.split(":")[0]
            tiling_reads.add(int(v))
            tiling_reads.add(int(w))
    return tiling_reads, tiling_path_data

if __name__ == "__main__":
    seqdb_prefix = sys.argv[1]
    tiling_path_fn = sys.argv[2]
    overlap_data_fn = sys.argv[3]
    total_chunk = int(sys.argv[4])
    my_chunk = int(sys.argv[5])

    seqdb = SequenceDatabase(
            f"{seqdb_prefix}.idx",
            f"{seqdb_prefix}.seqdb")


    # tiling_reads = set()
    tiling_reads, tiling_path_data = get_tiling_reads(tiling_path_fn, total_chunk, my_chunk)
    
    ovlps = get_ovlps(overlap_data_fn, tiling_reads)
    # ovlps = {}

    for ctg in tiling_path_data:
        chunks = []
        for i in range(0, len(tiling_path_data[ctg]), 32):
            chunk = tiling_path_data[ctg][i:i+32]
            chunks.append(chunk)

        template_cns_segments = []
      
        for chunk in chunks:
	    #template_seq, support_seq_ids = stiching_reads(chunk, seqdb)  
            template_seq, supporting_seq_ids = stiching_reads(chunk, seqdb, ovlps)
            print(ctg, "T:", len(template_seq), len(supporting_seq_ids), file=sys.stderr, flush=True)

            c_null = shimmer_ffi.NULL
            w = 128
            k = 16
            template_mmers = shimmer_ffi.new("mm128_v *")

            print("sketch start", file=sys.stderr, flush=True)
            shimmer4py.mm_sketch(c_null, template_seq, len(template_seq), w, k, 0, 0, template_mmers)
            print("sketch end", file=sys.stderr, flush=True)
            supporting_reads = []
            for r in supporting_seq_ids:
                rid2, tiling_read_direction, direction = r 
                if tiling_read_direction == 1:
                    s = seqdb.index_data[rid2].offset
                    l = seqdb.index_data[rid2].length
                    read_seq = shimmer_ffi.new("char[{}]".format(l))
                    shimmer4py.decode_biseq(seqdb.seqdb[s:s+l], read_seq, l, 1-direction) 
                else:
                    s = seqdb.index_data[rid2].offset
                    l = seqdb.index_data[rid2].length
                    read_seq = shimmer_ffi.new("char[{}]".format(l))
                    shimmer4py.decode_biseq(seqdb.seqdb[s:s+l], read_seq, l, direction) 
                offset = get_offset(template_mmers, read_seq, l, w, k)
                # print(rid2, offset)
                supporting_reads.append( (rid2, offset, read_seq) )

            print("cns start", file=sys.stderr, flush=True)
            cns = get_cns(template_seq, supporting_reads)
            print("cns end", file=sys.stderr, flush=True)
            template_cns_segments.append(cns)
            for _1, _2, read_seq in supporting_reads:
                shimmer_ffi.release(read_seq)
            
            #print(cns.decode("ascii"))
        print(ctg, "S:",len(template_cns_segments), file=sys.stderr, flush=True)

        template_cns_seqs  = []
        left_start = 0
        right_read = None
        if len(template_cns_segments) == 1:
            print(">"+ctg)
            print(template_cns_segments[0].decode("ascii"))
            continue
        for i in range(len(template_cns_segments)-1):
            left_read = template_cns_segments[i][left_start:]
            right_read = template_cns_segments[i+1]
            aln = get_shimmer_alns_from_seqs(left_read[-20000:], right_read[:20000], parameters={"w":80, "max_repeat":1}) # still dangeous in some repeat case
            aln = [_ for _ in aln if len(_[0]) > 5]
            aln.sort(key = lambda _: -len(_[0]))

            print("number of alns:", len(aln), file=sys.stderr, flush=True)
            if len(aln) == 1:
                smer0, smer1 = aln[0][0][0] # first mmer-pair
                left_end = smer0.pos_end + len(left_read) - 20000
                right_start = smer1.pos_end
            else:
                aln = get_shimmer_alns_from_seqs(left_read[-20000:], right_read[:20000], parameters={"w":32, "max_repeat":12})
                aln.sort(key = lambda _: -len(_[0]))
                if len(aln) > 0:
                    for j in range(len(aln)):
                        smer0, smer1 = aln[j][0][0] # first mmer-pair
                        left_end = smer0.pos_end + len(left_read) - 20000
                        right_start = smer1.pos_end
                        if left_end > right_start:
                            break
                    print("DEBUG_1_left", left_end, left_read.decode("ascii"), file = sys.stderr, flush=True)
                    print("DEBUG_1_right", right_start, right_read.decode("ascii"), file = sys.stderr, flush=True)
                else:
                    left_end = len(left_read)
                    right_start = 0
                    print("DEBUG_2_left", left_end, left_read, file=sys.stderr, flush=True)
                    print("DEBUG_2_right", right_start, right_read, file=sys.stderr, flush=Ture)
            template_cns_seqs.append(left_read[:left_end].decode("ascii"))
            left_start = right_start
        template_cns_seqs.append(right_read[right_start:].decode("ascii"))

        print(">"+ctg)
        print("".join(template_cns_seqs))
