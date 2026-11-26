set -euo pipefail
cd /staging/tmp/apple_minisat || exit 1
# loci: id fam chr start end len status dist_end where chrlen ; centromere: chr name start end
awk -F'\t' 'FNR==NR{cs[$1]=$3; ce[$1]=$4; next} {c=$3; mid=($4+$5)/2; if(mid>=cs[c]&&mid<=ce[c])d=0; else d=(mid<cs[c])?cs[c]-mid:mid-ce[c]; print $1"\t"$2"\t"c"\t"mid"\t"d"\t"$10}' cfg/centromere_hap1.tsv final_loci.tsv > loci_vs_centromere.tsv
echo "family n | inside the centromere range | within 1 Mb of it | uniform expectation (inside; within 1 Mb) | median distance (Mb)"
for g in $(cut -f2 final_loci.tsv | sort -u); do
 awk -F'\t' -v g=$g 'FNR==NR{cl[$1]=$4-$3; next} $2==g{n++; if($5==0)i++; if($5<=1000000)w++; ex_i+=cl[$3]/$6; ex_w+=(cl[$3]+2000000)/$6; d[n]=$5} END{asort(d); printf "%s n=%d | %d (%.0f%%) | %d (%.0f%%) | %.1f%%; %.1f%% | %.1f\n",g,n,i,100*i/n,w,100*w/n,100*ex_i/n,100*ex_w/n,d[int((n+1)/2)]/1e6}' <(awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$4}' cfg/centromere_hap1.tsv) loci_vs_centromere.tsv
done
echo "-- loci inside a centromere range (family chr locus-midpoint centromere)"
awk -F'\t' 'FNR==NR{cs[$1]=$3; ce[$1]=$4; nm[$1]=$2; next} $5==0{printf "%s %s %s %.0f in %d-%d\n",$2,nm[$3],$3,$4,cs[$3],ce[$3]}' cfg/centromere_hap1.tsv loci_vs_centromere.tsv | sort
