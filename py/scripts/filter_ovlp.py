import glob
reads = set()
with open("filtered_out_reads") as f:
    for r in f:
        r = r.strip()
        reads.add(r)

phasable_reads = set()
fn = "preads-ann.ovl"
with open(fn) as f:
    for r in f:
        r = r.strip().split()
        if r[0] in reads:
            continue
        if r[1] in reads:
            continue
        d = int(r[13])
        n = int(r[14])
        if d < 4 or 1.0*n/d >= 0.75:
            print(" ".join(r))
            phasable_reads.add(r[0])
            phasable_reads.add(r[1])

with open("phasable_reads", "w") as f:
    for rid in phasable_reads:
        print(rid, file=f)

unphased_reads = set()

fn ="preads-ann.ovl"
with open(fn) as f:
    for r in f:
        r = r.strip().split()
        if r[0] in reads:
            continue
        if r[1] in reads:
            continue
        if r[0] in phasable_reads:
            continue
        if r[1] in phasable_reads:
            continue
        unphased_reads.add(r[0])
        unphased_reads.add(r[1])
        #print(" ".join(r))
print("-")
with open("unphased_reads", "w") as f:
    for rid in unphased_reads:
        print(rid, file=f)
