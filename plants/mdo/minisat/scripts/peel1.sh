set -euo pipefail
export PATH=/staging/tmp/apple_minisat/env/bin:$PATH
mkdir -p /staging/tmp/apple_minisat/peel
cd /staging/tmp/apple_minisat/peel || exit 1
IN=/staging/tmp/apple_minisat/run_full/motifs/mdo.fa.clean.fas
# single-line FASTA with short ids m001..; map table
awk '/^>/{if(s!="")print id"\t"h"\t"s; n++; id=sprintf("m%03d",n); h=substr($0,2); s=""; next}{s=s $0} END{print id"\t"h"\t"s}' $IN > motifs.tsv
awk -F'\t' '{print ">"$1"\n"$3}' motifs.tsv > subj.fa
awk -F'\t' '{print ">"$1"\n"$3$3}' motifs.tsv > query2x.fa
makeblastdb -in subj.fa -dbtype nucl -out subj >/dev/null
blastn -task blastn -query query2x.fa -db subj -strand both -dust no -word_size 7 -evalue 1e-3 -max_target_seqs 1000 -max_hsps 3 \
  -outfmt '6 qseqid sseqid pident length qlen slen qstart qend sstart send sstrand' -num_threads 16 > allvsall.tsv
wc -l motifs.tsv allvsall.tsv
