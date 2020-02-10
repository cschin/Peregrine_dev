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

    cns = falcon4py.get_cns_from_align_tags(tags,
                                            aln_count,
                                            len(seq0),
                                            min_cov)

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
        return None
    for idx in range(aln.n):
        if aln.a[i].idx0.n > max_chain_n:
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

if __name__ == "__main__":
    seqdb_prefix = sys.argv[1]
    tiling_path_fn = sys.argv[2]
    overlap_data_fn = sys.argv[3]
    total_chunk = int(sys.argv[4])
    my_chunk = int(sys.argv[5])

    seqdb = SequenceDatabase(
            f"{seqdb_prefix}.idx",
            f"{seqdb_prefix}.seqdb")

    read_idx = {}
    with open("{}.idx".format(seqdb_prefix)) as f:
        for row in f:
            row = row.strip().split()
            rid, rname, rlen, offset = row
            rid = int(rid)
            rlen = int(rlen)
            offset = int(offset)
            read_idx.setdefault(rid, {})
            read_idx[rid]["name"] = rname
            read_idx[rid]["length"] = rlen
            read_idx[rid]["offset"] = offset

    tiling_path_data = {}
    tiling_reads = set()
    with open(tiling_path_fn) as f:
        for row in f:
            row = row.strip().split()
            ctg = row[0]
            ctg_hash = int(ctg[:-1]) 
            if ctg_hash % total_chunk != my_chunk % total_chunk:
                continue
            tiling_path_data.setdefault(row[0], [])
            tiling_path_data[row[0]].append(row)
            ctg_id, v, w, r, s, e, olen, idt, _1, _2 = row
            v = v.split(":")[0]
            w = w.split(":")[0]
            tiling_reads.add(int(v))
            tiling_reads.add(int(w))

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

    for ctg in tiling_path_data:
        chunks = []
        for i in range(0, len(tiling_path_data[ctg]), 24):
            chunk = tiling_path_data[ctg][i:i+24]
            chunks.append(chunk)

        template_cns_segments = []
        for chunk in chunks:
            template_seqs = []
            support_seq_ids = set()
            for row in chunk:
                ctg_id, v, w, r, s, e, olen, idt, _1, _2 = row
                r = int(r)
                tiling_read_direction = 0 if int(e) > int(s) else 1
                seq = seqdb.get_subseq_by_rid(r, direction=tiling_read_direction)
                template_seqs.append(seq)
                for rid2 in ovlps[r]:
                    ref_s, ref_e, s, e, direction = ovlps[r][rid2]
                    support_seq_ids.add( (rid2, tiling_read_direction, direction) ) 
            left_start = 0
            template_segments  = []
            for i in range(len(template_seqs)-1):
                left_read = template_seqs[i][left_start:]
                right_read = template_seqs[i+1]
                aln = get_shimmer_alns_from_seqs(left_read, right_read[:2500], parameters={"w":24, "max_repeat":12})
                aln.sort(key = lambda _: -len(_[0]))
                if len(aln) > 0:
                    smer0, smer1 = aln[0][0][0] # first mmer-pair
                    # print(i, len(aln[0][0]),
                    #       smer0.pos_end, smer0.direction, smer1.pos_end, smer1.direction, 
                    #       left_read[smer0.pos_end-16:smer0.pos_end],
                    #       right_read[smer1.pos_end-16:smer1.pos_end])
                    left_end = smer0.pos_end
                    right_start = smer1.pos_end
                else:
                    left_end = len(left_read)
                    right_start = 0
                    print("DEBUG_0_left", left_read)
                    print("DEBUG_0_right", right_read)
                template_segments.append(left_read[:left_end])
                left_start = right_start
            template_segments.append(right_read[right_start:])
            template_seq = b"".join(template_segments)

            c_null = shimmer_ffi.NULL
            w = 128
            k = 16
            template_mmers = shimmer_ffi.new("mm128_v *")
            shimmer4py.mm_sketch(c_null, template_seq, len(template_seq), w, k, 0, 0, template_mmers)
            supporting_reads = []
            for r in support_seq_ids:
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

            cns = get_cns(template_seq, supporting_reads)
            template_cns_segments.append(cns)
            for _1, _2, read_seq in supporting_reads:
                shimmer_ffi.release(read_seq)
            
            #print(cns.decode("ascii"))

        template_cns_seqs  = []
        left_start = 0
        for i in range(len(template_cns_segments)-1):
            left_read = template_cns_segments[i][left_start:]
            right_read = template_cns_segments[i+1]
            aln = get_shimmer_alns_from_seqs(left_read, right_read[:2500], parameters={"w":24, "max_repeat":12}) # still dangeous in some repeat case
            aln.sort(key = lambda _: -len(_[0]))
            
            # for a in aln:
            #     smer0, smer1 = a[0][0] # first mmer-pair
            #     print(len(aln), len(a[0]), smer0.pos_end, smer1.pos_end)
            # print()

            #smer0, smer1 = aln[0][0][0] # first mmer-pair
            if len(aln) > 0:
                smer0, smer1 = aln[0][0][0] # first mmer-pair
                # print(i, len(aln[0][0]),
                #       smer0.pos_end, smer0.direction, smer1.pos_end, smer1.direction, 
                #       left_read[smer0.pos_end-16:smer0.pos_end],
                #       right_read[smer1.pos_end-16:smer1.pos_end])
                left_end = smer0.pos_end
                right_start = smer1.pos_end
            else:
                left_end = len(left_read)
                right_start = 0
                print("DEBUG_1_left", left_read)
                print("DEBUG_1_right", right_read)
            # print(i, len(aln[0][0]),
            #       smer0.pos_end, smer0.direction, smer1.pos_end, smer1.direction, 
            #       left_read[smer0.pos_end-16:smer0.pos_end],
            #       right_read[smer1.pos_end-16:smer1.pos_end])
            left_end = smer0.pos_end
            right_start = smer1.pos_end
            template_cns_seqs.append(left_read[:left_end].decode("ascii"))
            left_start = right_start

        template_cns_seqs.append(right_read[right_start:].decode("ascii"))
        print(">"+ctg)
        print("".join(template_cns_seqs))
