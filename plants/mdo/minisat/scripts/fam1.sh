set -euo pipefail
cd /staging/tmp/apple_minisat || exit 1
rm -rf repo2 && tar xzf repo2.tgz
mkdir -p cfg/fam
# one single-record FASTA per peeled group, header = family name
awk '/^>cons_/{split(substr($0,7),a,"_"); g=a[1]; getline s; L=length(s); print ">"g"_"L"\n"s > "cfg/fam/"g"_"L".fa"; printf "%s_%d\t%d\n", g, L, L}' peel/consensus_all.fa | sort > cfg/fam/families.tsv
cat cfg/fam/families.tsv
{
echo "# Apple de novo minisatellite scan, pass 2: families = the seven groups peeled from the stage-04 motifs (peel/ on DRAGEN)."
echo "run_name: apple_minisat_fam"
echo "out_dir: /staging/tmp/apple_minisat/run_fam"
echo "genomes:"
echo "  mdo: /staging/tmp/apple_sines/genome/apple.fna"
echo "families:"
while IFS=$'\t' read -r name mono; do echo "  $name: {consensus: /staging/tmp/apple_minisat/cfg/fam/$name.fa, mono: $mono}"; done < cfg/fam/families.tsv
echo "anchor:"
echo "  independent: true"
echo "family_assignment: /staging/tmp/apple_minisat/cfg/cluster_family_fam.tsv"
} > cfg/apple_minisat_fam.yaml
: > cfg/cluster_family_fam.tsv
cat cfg/apple_minisat_fam.yaml
export PATH=/staging/tmp/apple_minisat/env/bin:/staging/tmp/apple_minisat/repo2/bin:$PATH
export TMPDIR=/staging/tmp/apple_minisat/tmp
satpipe run /staging/tmp/apple_minisat/cfg/apple_minisat_fam.yaml --to 08_cross > fam_to08.log 2>&1 || true
grep -v '^\[OK' fam_to08.log | tail -30 | cut -c1-200
