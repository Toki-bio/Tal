set -euo pipefail
cd /staging/tmp/apple_minisat || exit 1
export CONDA_PKGS_DIRS=/staging/tmp/apple_minisat/pkgs
export TMPDIR=/staging/tmp/apple_minisat/tmp
mkdir -p "$TMPDIR"
/staging/miniconda3/bin/mamba env create -y -p /staging/tmp/apple_minisat/env -f /staging/tmp/apple_minisat/environment.yml
touch /staging/tmp/apple_minisat/env.DONE
