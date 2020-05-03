#for i in `seq 1 20`; do echo "./shmr_overlap_annotate -p ../0-seqdb/seq_dataset -l ./preads.ovl -t 20 -c $i > preads-ann.ovl.$i"; done | parallel -j 20
#cat preads-ann.ovl.* | awk '$NF == 1 {print $1}' | sort -u > filtered_out_reads
python filter_ovlp.py > tmp2.ovl
/usr/bin/time pypy3 ./ovlp_to_graph_exp.py --overlap-file tmp2.ovl --min_len 8000 --min_idt 96 >& asm.log
/usr/bin/time pypy3 ./graph_to_path.py >& to_path.log
for i in `seq 1 20`; do echo "/usr/bin/time python tiling_path_ec3.py ../0-seqdb/seq_dataset p_ctg_tiling_path tmp2.ovl 20 $i  > tmp_p_dbg_chunk_$i.fa 2> log.p.$i"; done | parallel -j 20 
cat tmp_p_dbg_chunk_*.fa > p_phase_ctg.fa
for i in `seq 1 20`; do echo "/usr/bin/time python tiling_path_ec3.py ../0-seqdb/seq_dataset a_ctg_tiling_path tmp2.ovl 20 $i  > tmp_a_dbg_chunk_$i.fa 2> log.a.$i"; done | parallel -j 20 
cat tmp_a_dbg_chunk_*.fa > a_phase_ctg.fa
cat p_phase_ctg.fa a_phase_ctg.fa > pa_phase_ctg.fa


