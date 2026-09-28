source ~/miniforge3/etc/profile.d/conda.sh; conda activate sinederella
O=~/tmp/recons; rm -rf $O; mkdir -p $O; cd $O
SECONDS=0
python3 ~/tmp/inst_new/recons.py $O 2> $O/recons.err; echo "collect: ${SECONDS}s"; tail -3 $O/recons.err
grep -c ">" recreated.fa
mafft --localpair --maxiterate 1000 --adjustdirectionaccurately --reorder --preservecase --thread 16 --quiet recreated.fa > recreated.aln.fa 2> mafft.err; echo "mafft exit=$? ${SECONDS}s"
grep -c "^>_R_" recreated.aln.fa
python3 ~/tmp/inst_new/nearest.py recreated.aln.fa members.tsv nearest.tsv
echo "--- rebuilt consensuses whose closest reference is NOT their own query:"
awk -F'\t' 'NR>1 && $8 != $4' nearest.tsv | cut -f1,5,8-13 | column -t
echo "--- counts by own family -> closest"; awk -F'\t' 'NR>1{print $4" -> "$8}' nearest.tsv | sort | uniq -c | sort -rn
echo "--- skipped"; cat skipped.tsv | cut -f1-4 | column -t | head -60
