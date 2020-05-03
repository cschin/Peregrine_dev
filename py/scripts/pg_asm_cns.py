#!/usr/bin/env python3

import mmap
import sys
import re
from collections import OrderedDict
from peregrine._falcon4py import ffi
from peregrine._falcon4py import lib as falcon
from peregrine._shimmer4py import lib as shimmer

low_coverage_m = re.compile(rb"[acgt]+")

## No option parsing at thie moment, perhaps letter

read_db_prefix = sys.argv[1]
ref_db_prefix = sys.argv[2]
read_to_contig_map = sys.argv[3]
total_chunks = int(sys.argv[4])
my_chunk = int(sys.argv[5])


f = open("{}.seqdb".format(read_db_prefix), "rb")
seqdb = mmap.mmap(f.fileno(), 0, flags=mmap.MAP_SHARED, prot=mmap.PROT_READ)

f = open("{}.seqdb".format(ref_db_prefix), "rb")
refdb = mmap.mmap(f.fileno(), 0, flags=mmap.MAP_SHARED, prot=mmap.PROT_READ)

read_idx = {}
with open("{}.idx".format(read_db_prefix)) as f:
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


ref_idx = {}
with open("{}.idx".format(ref_db_prefix)) as f:
    for row in f:
        row = row.strip().split()
        rid, rname, rlen, offset = row
        rid = int(rid)
        rlen = int(rlen)
        offset = int(offset)
        ref_idx.setdefault(rid, {})
        ref_idx[rid]["name"] = rname
        ref_idx[rid]["length"] = rlen
        ref_idx[rid]["offset"] = offset

contig_to_read_map = OrderedDict()
with open(read_to_contig_map) as f:
    for row in f:
        row = row.strip().split()
        row = tuple(int(c) for c in row)
        ctg_id = row[0]
        if (my_chunk % total_chunks) != (ctg_id % total_chunks):
            continue
        contig_to_read_map.setdefault(ctg_id, [])
        contig_to_read_map[ctg_id].append(row)

rng = ffi.new("aln_range[1]")


# TODO: we need to refactor this loop
for ctg in contig_to_read_map:
    print("-\n", "ctg {}".format(ref_idx[ctg]["name"]), file=sys.stderr)
    contig_to_read_map[ctg].sort(key=lambda x: x[1])
    read_map_groups = []
    left_anchor = 1000
    map_group = []

    for row in contig_to_read_map[ctg]:
        ref_p1 = row[1]
        if ref_p1 - left_anchor < 50000:
            map_group.append(row)
        else:
            if ref_p1 - left_anchor < 100000:
                read_map_groups.append([left_anchor, ref_p1, map_group])
            else:
                read_map_groups.append([left_anchor, ref_p1, []])
            map_group = []
            left_anchor = ref_p1

    if ref_idx[ctg]["length"] - left_anchor < 100000:  #current max template size for consensus
        if ref_idx[ctg]["length"] - left_anchor > 1000:
            read_map_groups.append((left_anchor,
                                    ref_idx[ctg]["length"],
                                    map_group))
        elif len(read_map_groups) > 0:
            read_map_groups[-1][1] = ref_idx[ctg]["length"]
            read_map_groups[-1][2].extend(map_group)
        else:
            read_map_groups.append((left_anchor, ref_idx[ctg]["length"], []))
    else:
        read_map_groups.append((left_anchor, ref_idx[ctg]["length"], []))

    print("ctg {}".format(ref_idx[ctg]["name"]),
          len(read_map_groups),
          file=sys.stderr)

    #if len(read_map_groups) <= 2: #ignore short contig for now
    #    continue

    cns_segments = []
    j = 0
    for left, right, mapped in read_map_groups:
        print(f"--\n sg{j:03d}", left, right, right-left, len(mapped), file=sys.stderr)

        j += 1
        left = left-1000
        assert(left >= 0)
        rmap = {}

        for d in mapped:
            #print(d)
            read_id = d[3]
            read_offset = d[1] - d[4]
            read_strand = d[6]
            rmap.setdefault((read_id, read_strand), [])
            rmap[(read_id, read_strand)].append(read_offset)

        reads = []

        for (read_id, read_strand), v in rmap.items():
            v.sort()
            v_current = v[0]
            reads.append((read_id, read_strand, v_current - left, len(v)))
            print( (read_id, read_strand), v_current, file=sys.stderr);
            for vv in v:
                if vv > v_current + 50:
                    v_current = vv
                    reads.append((read_id, read_strand, v_current - left, len(v)))
                    print( (read_id, read_strand), v_current, file=sys.stderr);


        reads.sort(key=lambda x: x[2])
        s = ref_idx[ctg]["offset"] + left
        ref_len = right-left

        bseq0 = refdb[s:s+ref_len]

        ref_seq = ffi.new("char[{}]".format(ref_len))

        shimmer.decode_biseq(bseq0, ref_seq, ref_len, 0)

        tags = ffi.new("align_tags_t * [{}]".format(len(reads)+1))

        # need a back bone for some boundary case
        aln = falcon.align(ref_seq, ref_len,
                           ref_seq, ref_len,
                           50, 1)
        rng[0].s1 = aln.aln_q_s
        rng[0].e1 = aln.aln_q_e
        rng[0].s2 = aln.aln_t_s
        rng[0].e2 = aln.aln_t_e
        tag = falcon.get_align_tags(aln.q_aln_str,
                                    aln.t_aln_str,
                                    aln.aln_str_size,
                                    rng, 0, 0)
        aln_count = 0
        tags[aln_count] = tag
        aln_count += 1
        falcon.free_alignment(aln)

        aln_base = 0
        for d in reads:
            #print(d)
            read_id = d[0]
            read_strand = d[1]
            read_shift = int(d[2])
            s = read_idx[read_id]["offset"]
            read_len = read_idx[read_id]["length"]
            bseq1 = seqdb[s:s+read_len]
            read_seq = ffi.new("char[{}]".format(read_len))
            shimmer.decode_biseq(bseq1, read_seq, read_len, read_strand)

            aligned = False
            t_offset = 0
            if read_shift < 0:
                aln = falcon.align(read_seq[abs(read_shift):read_len],
                                   read_len - abs(read_shift),
                                   ref_seq,
                                   ref_len,
                                   150, 1)

                if abs(abs(aln.aln_q_e-aln.aln_q_s) -
                       (read_len - abs(read_shift))) < 48:
                    aligned = True

                    rng[0].s1 = aln.aln_q_s
                    rng[0].e1 = aln.aln_q_e
                    rng[0].s2 = aln.aln_t_s
                    rng[0].e2 = aln.aln_t_e
                    t_offset = 0
                else:
                    falcon.free_alignment(aln)
            else:
                aln = falcon.align(read_seq,
                                   read_len,
                                   ref_seq[read_shift:ref_len],
                                   ref_len-read_shift,
                                   150, 1)

                if abs(abs(aln.aln_q_e-aln.aln_q_s)-read_len) < 48 or \
                   abs(ref_len-read_shift-abs(aln.aln_q_e-aln.aln_q_s)) < 48:
                    aligned = True
                    rng[0].s1 = aln.aln_q_s
                    rng[0].e1 = aln.aln_q_e
                    rng[0].s2 = aln.aln_t_s
                    rng[0].e2 = aln.aln_t_e
                    t_offset = read_shift
                else:
                    falcon.free_alignment(aln)
            if aligned:
                print(f"{read_id} is algined",
                      rng[0].s1 , rng[0].e1, rng[0].s2, rng[0].e2, file=sys.stderr)
                # print(ffi.string(aln.q_aln_str), file=sys.stderr)
                # rint(ffi.string(aln.t_aln_str), file=sys.stderr)
                sys.stderr.flush()
                tag = falcon.get_align_tags(aln.q_aln_str,
                                            aln.t_aln_str,
                                            aln.aln_str_size,
                                            rng, 0, t_offset)
                tags[aln_count] = tag
                aln_count += 1
                aln_base += abs(rng[0].e2 - rng[0].s2)
                falcon.free_alignment(aln)
            ffi.release(read_seq)
        aln_cov = aln_base/ref_len
        print(f"aln_count:{aln_count}, aln_base: {aln_base}, aln_cov: {aln_cov}", file=sys.stderr)

        if aln_base/ref_len < 5:
            cns_seq = ffi.string(ref_seq)
            cns_seq = cns_seq.lower()
        else:
            cns = falcon.get_cns_from_align_tags(tags,
                                                 aln_count, len(ref_seq), 1)
            cns_seq = ffi.string(cns.sequence)
            falcon.free_consensus_data(cns)

        cns_segments.append(cns_seq)

        for i in range(aln_count):
            falcon.free_align_tags(tags[i])
        ffi.release(tags)
        ffi.release(ref_seq)

    s0 = cns_segments[0]
    stiched_segments = [s0]
    p = 0
    for s1 in cns_segments[1:]:
        aln = falcon.align(s0[-1000:], 1000,
                           s1[:1050], 1050, 400, 0)
        # print(aln.aln_q_s, aln.aln_q_e, aln.aln_t_s, aln.aln_t_e, aln.dist)
        if aln.aln_q_e < 1000:
            stiched_segments[-1] = stiched_segments[-1][:-(1000-aln.aln_q_e)]

        stiched_segments.append(s1[aln.aln_t_e:])
        p += len(s1[aln.aln_t_e:])
        print("stiching point:", p, file=sys.stderr)
        print("aln.aln_q_e:", aln.aln_q_e, file=sys.stderr)
        print("aln.aln_t_e:", aln.aln_t_e, file=sys.stderr)
        # print(ffi.string(aln.q_aln_str), file=sys.stderr)
        # print(ffi.string(aln.t_aln_str), file=sys.stderr)
        s0 = s1
        falcon.free_alignment(aln)

    contig = b"".join(stiched_segments)
    low_cov_count = 0
    for m in low_coverage_m.finditer(contig):
        s_pos, e_pos = m.span()
        low_cov_count += e_pos - s_pos

    # we might want to collect low coverage consensus in a different file
    # in the future
    if 1.0 * low_cov_count/len(contig) < 0.5:
        print(">{}".format(ref_idx[ctg]["name"]))
        print(contig.decode("ascii"))

ffi.release(rng)
