#!/usr/bin/env bash
# after_pull.sh [code ...] - run after pulling republished bat reports into Tal.
# A republish (SINEderella publish_run.sh) regenerates <code>/report.html from scratch, which drops the
# blocks Tal adds: back navigation + species line (build_lineages.add_back_nav) and the SINE Tribes
# section (inject_tribes.py). Both are idempotent. Happened 2026-09-28: all 25 pages lost them.
# Does NOT rebuild chiroptera.html (build_lineages.py main) - that page has hand-added cards.
cd "$(dirname "$0")/.." || exit 1
codes="${*:-$(awk -F'\t' 'NR>1 && $1!="-"{print $1}' chiroptera/genomes.tsv)}"
PYTHONIOENCODING=utf-8 python3 - $codes <<'PY'
import os, sys
sys.path.insert(0, "chiroptera")
import build_lineages as b
n = sum(b.add_back_nav(c + "/report.html") for c in sys.argv[1:] if os.path.exists(c + "/report.html"))
print(n, "reports: navigation + species line")
PY
python3 chiroptera/inject_tribes.py $codes | tail -1
