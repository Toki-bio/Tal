cd /staging/tmp/apple_minisat || exit 1
echo "chromosome | n G02 loci | nearest distance to centromere range (kb) | centromere range (Mb) | chromosome length (Mb)"
awk -F'\t' 'FNR==NR{cs[$1]=$3; ce[$1]=$4; nm[$1]=$2; next} $2=="G02_170"{n[$3]++; if(!($3 in m)||$5<m[$3])m[$3]=$5; L[$3]=$6} END{for(c in n)printf "%s %s | %d | %.0f | %.2f-%.2f | %.1f\n",nm[c],c,n[c],m[c]/1000,cs[c]/1e6,ce[c]/1e6,L[c]/1e6}' cfg/centromere_hap1.tsv loci_vs_centromere.tsv | sort
echo "chromosomes with a centromere range: 17; G02 on 8 chromosomes"
