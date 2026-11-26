set -euo pipefail
export PATH=/staging/tmp/apple_minisat/env/bin:$PATH
cd /staging/tmp/apple_minisat/peel || exit 1
A=/staging/tmp/apple_sines
cat $A/sd_launch9.peel7.fa $A/repbase_plant84.fa > sinebank.fa
grep -c '>' sinebank.fa
makeblastdb -in sinebank.fa -dbtype nucl -out sinebank >/dev/null
awk '/^>/{if(h!="")print ">"h"\n"s s; h=substr($0,2); s=""; next}{s=s $0} END{print ">"h"\n"s s}' consensus_all.fa > cons2x.fa
echo "== consensus (doubled) vs SINE bank (7 owner peeled consensuses + 84 Repbase land-plant SINEs), best hit per pair, both strands"
blastn -task blastn -query cons2x.fa -db sinebank -strand both -dust no -word_size 7 -evalue 1e-4 -outfmt '6 qseqid sseqid pident length qlen slen sstrand' | sort -k1,1 -k4,4nr | awk '{c[$1]++; if(c[$1]<=3)printf "%s -> %s pid=%.1f aln=%d (subject len %d) %s\n",$1,$2,$3,$4,$6,$7}'
echo "== SINE_discriminator-style: any tRNA-like? G03 T-arm"
awk '/^>cons_G03/{getline s; print s}' consensus_all.fa | awk '{print index($0,"GGTTCGA"), index($0,"TTCGACTCC")}'
echo "== G04 translation, 3 forward + 3 reverse frames (longest ORF-ish stretch without stop)"
awk '/^>cons_G04/{getline s; print ">G04\n"s s}' consensus_all.fa > g04.fa
transeq -sequence g04.fa -frame 6 -clean -outseq g04.pep >/dev/null 2>&1 || transeq -sequence g04.fa -frame 6 -outseq g04.pep >/dev/null 2>&1
awk '/^>/{h=$0; next}{n=split($0,a,"*"); mx=0; for(i=1;i<=n;i++)if(length(a[i])>mx){mx=length(a[i]);b=a[i]} print h, "stops:",n-1, "longest stop-free:",mx; print b}' g04.pep
echo "== files in the SINE run results (for locus overlap)"
ls $A/sd_peel9/run_excl_20261126_031836/results | head -20
