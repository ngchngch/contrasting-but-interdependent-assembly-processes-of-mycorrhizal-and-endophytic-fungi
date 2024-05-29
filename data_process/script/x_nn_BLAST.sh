
_T=n
str='XXXXXXX===nt database path===XXXXXXX'
str="$(eval echo $str | cut -d '.' -f1)"

fas='05_Merge_sequence_output/Plant/OTUseq_0.97.fasta'
out='06_Taxonomy_annotaion/OTUseq_0.97'
nn=5
iden=95
cov=20
thr=48


./bin/export_BLAST_query $fas $out $nn

case "${_T}" in
	[N/n])
		### hsp is the sequential base that have high homology
		echo ""
		echo ""
		echo "./bin/ncbi-blast-2.13.0+/bin/blastn -db $str -query $fas -out $out -evalue 10 -perc_identity $iden -qcov_hsp_perc $cov -outfmt '6 qseqid sseqid evalue score length pident sstrand qcovs' -max_target_seqs $nn -num_threads $thr -strand 'both'"
		echo ""
		echo ""
		./bin/ncbi-blast-2.13.0+/bin/blastn -db $str -query $fas -out $out -evalue 10 -perc_identity $iden -qcov_hsp_perc $cov -outfmt '6 qseqid sseqid evalue score length pident sstrand qcovs' -max_target_seqs $nn -num_threads $thr -strand 'both'
		
		./bin/x_nn $out $nn
		
		echo ""
		echo "******************************"
		echo "          Completed!"
		echo "******************************"
		echo ""
	;;
	[P/p])
		### hsp is the sequential base that have high homology
		echo ""
		echo ""
		echo "./bin/ncbi-blast-2.13.0+/bin/blastp -db $str -query $fas -out $out -evalue 10 -perc_identity $iden -qcov_hsp_perc $cov -outfmt '6 qseqid sseqid evalue score length pident sstrand qcovs' -max_target_seqs $nn -num_threads $thr -strand 'both'"
		echo ""
		echo ""
		./bin/ncbi-blast-2.13.0+/bin/blastn -db $str -query $fas -out $out -evalue 10 -perc_identity $iden -qcov_hsp_perc $cov -outfmt '6 qseqid sseqid evalue score length pident sstrand qcovs' -max_target_seqs $nn -num_threads $thr -strand 'both'
		
		./bin/x_nn $out $nn
		
		echo ""
		echo "******************************"
		echo "          Completed!"
		echo "******************************"
		echo ""
	;;
	[X/x])
		### hsp is the sequential base that have high homology
		echo ""
		echo ""
		echo "./bin/ncbi-blast-2.13.0+/bin/blastx -db $str -query $fas -out $out -evalue 10 -perc_identity $iden -qcov_hsp_perc $cov -outfmt '6 qseqid sseqid evalue score length pident sstrand qcovs' -max_target_seqs $nn -num_threads $thr -strand 'both'"
		echo ""
		echo ""
		./bin/ncbi-blast-2.13.0+/bin/blastn -db $str -query $fas -out $out -evalue 10 -perc_identity $iden -qcov_hsp_perc $cov -outfmt '6 qseqid sseqid evalue score length pident sstrand qcovs' -max_target_seqs $nn -num_threads $thr -strand 'both'
		
		./bin/x_nn $out $nn
		
		echo ""
		echo "******************************"
		echo "          Completed!"
		echo "******************************"
		echo ""
	;;esac


		
