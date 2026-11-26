set -euo pipefail
cd /staging/tmp/apple_minisat || exit 1
awk -F'\t' '$2!="-" && $2!=""{printf "%s\t%s\tpeeled group majority (%s of %s members)\n",$1,$2,$3,$4}' cfg/cluster_map_draft.tsv > cfg/cluster_family_fam.tsv
cat cfg/cluster_family_fam.tsv
export PATH=/staging/tmp/apple_minisat/env/bin:/staging/tmp/apple_minisat/repo2/bin:$PATH
export TMPDIR=/staging/tmp/apple_minisat/tmp
satpipe run /staging/tmp/apple_minisat/cfg/apple_minisat_fam.yaml > fam_full.log 2>&1
touch /staging/tmp/apple_minisat/fam_full.DONE
