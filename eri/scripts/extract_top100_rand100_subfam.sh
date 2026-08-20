#!/bin/bash
# extract_top100_rand100_subfam.sh
#
# Per-subfamily result extraction for the Tal SINE site, generalized from
# _sinederella/patch/run_subfam_per_sf.sh (which only covers the "subfam"
# variant and is hardcoded to the saq/ccr runs). This script produces all
# three variants used on the existing species pages (saq/ccr/toc/teu/gpy/dmo):
#
#   {species}_{sf}_top100.aln.fa   — 100 highest-bitscore copies + consensus, mafft-aligned
#   {species}_{sf}_rand100.aln.fa  — 100 randomly-sampled copies + consensus, mafft-aligned
#   {species}_{sf}_subfam.aln.fa   — SubFam re-clustering of up to N_SAMPLE copies
#                                     (chunk-consensuses + sf-consensus), same
#                                     method/parameters as run_subfam_per_sf.sh
#
# DELIBERATELY KEPT SEPARATE from _sinederella/patch/ (shared, used by other
# species) — this copy is scoped to eri/ only, per explicit instruction not to
# let eri-specific work spread into other Tal species pages. If this script
# proves useful beyond eri, promoting it into _sinederella/patch/ (as a
# parameterized replacement for run_subfam_per_sf.sh) is a separate decision,
# not made here.
#
# Design choices made without an existing example to copy exactly (WORKFLOW.md
# and step6_report.py reference top100/rand100 files but no generation script
# for them was found in this repo) — confirm these before trusting the output
# if anything looks off compared to the other species' existing files:
#   - top100 ranking: bitscore parsed from assigned.fasta headers
#     (>seqID|subfamily|bitscore), matching SINEderella's own primary metric.
#   - Both top100 and rand100 use the same mafft parameters as the "subfam"
#     variant's final alignment step (--localpair --maxiterate 1000 --ep 0.123
#     --nuc --reorder --preservecase), for visual/structural consistency
#     across the three alignment types on a report page.
#   - rand100 uses the same fixed seed (42) as run_subfam_per_sf.sh's sampling,
#     for reproducibility.
#
# Usage:
#   extract_top100_rand100_subfam.sh <RUN_ROOT> <SPECIES_CODE> <OUT_DIR>
#
# <RUN_ROOT>      SINEderella run directory (contains consensuses.clean.fa,
#                  step2/step2_output/assigned.fasta)
# <SPECIES_CODE>  short code used in output filenames, e.g. "eri"
# <OUT_DIR>       where to write the *.aln.fa files (created if missing)

set -euo pipefail

RUN_ROOT="${1:?usage: $0 <RUN_ROOT> <SPECIES_CODE> <OUT_DIR>}"
SPECIES="${2:?usage: $0 <RUN_ROOT> <SPECIES_CODE> <OUT_DIR>}"
OUT_DIR="${3:?usage: $0 <RUN_ROOT> <SPECIES_CODE> <OUT_DIR>}"

ASSIGNED="$RUN_ROOT/step2/step2_output/assigned.fasta"
CONS_FA="$RUN_ROOT/consensuses.clean.fa"
SUBFAM_BIN="/data/V/toki/bin/SINEderella/SubFam"
BINSIZE=50
N_SAMPLE=10000       # subfam variant: max copies fed to SubFam clustering
MIN_COPIES_SUBFAM=400  # subfam variant: skip subfamilies smaller than this
TOPN=100
RANDN=100
SEED=42
THREADS="$(nproc)"

[[ -s "$ASSIGNED" ]] || { echo "ERROR: missing/empty $ASSIGNED" >&2; exit 1; }
[[ -s "$CONS_FA"  ]] || { echo "ERROR: missing/empty $CONS_FA" >&2; exit 1; }

# Portable python3 selection (matches the fix already applied to
# step4_plots.sh — don't hardcode a single interpreter path that may not
# exist or lack needed stdlib features on every system this runs on).
PYTHON3="python3"
command -v "$PYTHON3" >/dev/null 2>&1 || {
  for cand in python3.13 python3.12 python3.11 python3.10 python3.9; do
    command -v "$cand" >/dev/null 2>&1 && { PYTHON3="$cand"; break; }
  done
}
command -v "$PYTHON3" >/dev/null 2>&1 || { echo "ERROR: no python3 found" >&2; exit 1; }

mkdir -p "$OUT_DIR"
WORK="$OUT_DIR/_work"
rm -rf "$WORK"
mkdir -p "$WORK"

echo "[$(date '+%H:%M:%S')] Splitting $ASSIGNED by subfamily, ranking by bitscore..."

# ---- Split assigned.fasta into per-subfamily FASTA, sorted by bitscore desc ----
# Header format: >seqID|subfamily|bitscore
"$PYTHON3" - "$ASSIGNED" "$WORK" <<'PYEOF'
import sys, os
assigned, work = sys.argv[1], sys.argv[2]

by_sf = {}
name = None
seq_lines = []

def flush():
    if name is None:
        return
    parts = name.split('|')
    if len(parts) < 3:
        return
    sf = parts[1]
    try:
        bs = int(parts[2])
    except ValueError:
        bs = 0
    by_sf.setdefault(sf, []).append((bs, name, ''.join(seq_lines)))

with open(assigned) as f:
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('>'):
            flush()
            name = line[1:]
            seq_lines = []
        else:
            seq_lines.append(line)
flush()

for sf, copies in sorted(by_sf.items()):
    copies.sort(key=lambda t: t[0], reverse=True)
    out = os.path.join(work, sf + '.sorted.fasta')
    with open(out, 'w') as fh:
        for bs, hdr, seq in copies:
            fh.write('>{}\n{}\n'.format(hdr, seq))
    print('{}: {} copies'.format(sf, len(copies)), flush=True)
PYEOF

# ---- Extract individual subfamily consensus sequences ----
"$PYTHON3" - "$CONS_FA" "$WORK" <<'PYEOF'
import sys
from pathlib import Path
cons_fa, work = sys.argv[1], sys.argv[2]
records = {}
name = None
seq_lines = []
with open(cons_fa) as f:
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if name:
                records[name] = ''.join(seq_lines)
            name = line[1:].split()[0]
            seq_lines = []
        else:
            seq_lines.append(line)
if name:
    records[name] = ''.join(seq_lines)
for sf_name, seq in records.items():
    out = Path(work) / (sf_name + '.cons.fasta')
    out.write_text('>{}\n{}\n'.format(sf_name + '_CONSENSUS', seq))
PYEOF

MAFFT_OPTS=(--thread "$THREADS" --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --preservecase --quiet)

for sorted_fa in "$WORK"/*.sorted.fasta; do
    [[ -f "$sorted_fa" ]] || continue
    sf="$(basename "${sorted_fa%.sorted.fasta}")"
    cons="$WORK/${sf}.cons.fasta"
    n_total="$(grep -c '^>' "$sorted_fa")"

    if [[ ! -f "$cons" ]]; then
        echo "[$(date '+%H:%M:%S')] WARN  $sf: no consensus fasta found, skipping"
        continue
    fi

    echo "[$(date '+%H:%M:%S')] $sf: $n_total total assigned copies"

    # ---- top100 ----
    head_n=$(( n_total < TOPN ? n_total : TOPN ))
    "$PYTHON3" - "$sorted_fa" "$head_n" "$WORK/${sf}.top${TOPN}.fasta" <<'PYEOF'
import sys
fa, n, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
recs = []
name = None
seq = []
with open(fa) as f:
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if name: recs.append((name, ''.join(seq)))
            name = line[1:]; seq = []
        else:
            seq.append(line)
if name: recs.append((name, ''.join(seq)))
with open(out, 'w') as fh:
    for nm, sq in recs[:n]:
        fh.write('>{}\n{}\n'.format(nm, sq))
PYEOF
    cat "$WORK/${sf}.top${TOPN}.fasta" "$cons" > "$WORK/${sf}.top${TOPN}.combined.fasta"
    mafft "${MAFFT_OPTS[@]}" "$WORK/${sf}.top${TOPN}.combined.fasta" \
        > "$OUT_DIR/${SPECIES}_${sf}_top100.aln.fa"
    echo "[$(date '+%H:%M:%S')]   OK top100 ($head_n copies + consensus)"

    # ---- rand100 (seeded) ----
    rand_n=$(( n_total < RANDN ? n_total : RANDN ))
    "$PYTHON3" - "$sorted_fa" "$rand_n" "$SEED" "$WORK/${sf}.rand${RANDN}.fasta" <<'PYEOF'
import sys, random
fa, n, seed, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
recs = []
name = None
seq = []
with open(fa) as f:
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if name: recs.append((name, ''.join(seq)))
            name = line[1:]; seq = []
        else:
            seq.append(line)
if name: recs.append((name, ''.join(seq)))
rng = random.Random(seed)
sampled = rng.sample(recs, n)
with open(out, 'w') as fh:
    for nm, sq in sampled:
        fh.write('>{}\n{}\n'.format(nm, sq))
PYEOF
    cat "$WORK/${sf}.rand${RANDN}.fasta" "$cons" > "$WORK/${sf}.rand${RANDN}.combined.fasta"
    mafft "${MAFFT_OPTS[@]}" "$WORK/${sf}.rand${RANDN}.combined.fasta" \
        > "$OUT_DIR/${SPECIES}_${sf}_rand100.aln.fa"
    echo "[$(date '+%H:%M:%S')]   OK rand100 ($rand_n copies + consensus)"

    # ---- subfam (SubFam re-clustering, same method as run_subfam_per_sf.sh) ----
    if (( n_total < MIN_COPIES_SUBFAM )); then
        echo "[$(date '+%H:%M:%S')]   SKIP subfam: $n_total < $MIN_COPIES_SUBFAM"
        continue
    fi
    sample_n=$(( n_total < N_SAMPLE ? n_total : N_SAMPLE ))
    "$PYTHON3" - "$sorted_fa" "$sample_n" "$SEED" "$WORK/${sf}.subfamsample.fasta" <<'PYEOF'
import sys, random
fa, n, seed, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
recs = []
name = None
seq = []
with open(fa) as f:
    for line in f:
        line = line.rstrip('\n')
        if line.startswith('>'):
            if name: recs.append((name, ''.join(seq)))
            name = line[1:]; seq = []
        else:
            seq.append(line)
if name: recs.append((name, ''.join(seq)))
rng = random.Random(seed)
sampled = rng.sample(recs, n)
with open(out, 'w') as fh:
    for nm, sq in sampled:
        fh.write('>{}\n{}\n'.format(nm, sq))
PYEOF
    sf_work="$WORK/_subfam_${sf}"
    rm -rf "$sf_work"; mkdir -p "$sf_work"
    cp "$WORK/${sf}.subfamsample.fasta" "$sf_work/input.fasta"
    ( cd "$sf_work" && "$SUBFAM_BIN" input.fasta "$BINSIZE" > subfam.log 2>&1 ) || {
        echo "[$(date '+%H:%M:%S')]   WARN subfam: SubFam failed, see $sf_work/subfam.log"
        continue
    }
    if [[ ! -f "$sf_work/input.clw" ]]; then
        echo "[$(date '+%H:%M:%S')]   WARN subfam: no input.clw produced (see $sf_work/subfam.log)"
        continue
    fi
    if head -1 "$sf_work/input.clw" | grep -q "^CLUSTAL"; then
        awk '
          /^CLUSTAL/ || /^$/ { next }
          NF >= 2 {
            name=$1; seq=$2
            gsub(/[-.]/, "", seq)
            if (!(name in s)) { ord[++n]=name }
            s[name] = s[name] seq
          }
          END { for (i=1;i<=n;i++) { print ">" ord[i]; print s[ord[i]] } }
        ' "$sf_work/input.clw" > "$sf_work/input_reps.fasta"
    else
        cp "$sf_work/input.clw" "$sf_work/input_reps.fasta"
    fi
    cat "$sf_work/input_reps.fasta" "$cons" > "$sf_work/combined.fasta"
    mafft "${MAFFT_OPTS[@]}" "$sf_work/combined.fasta" \
        > "$OUT_DIR/${SPECIES}_${sf}_subfam.aln.fa"
    echo "[$(date '+%H:%M:%S')]   OK subfam ($sample_n copies clustered + consensus)"
    rm -rf "$sf_work"
done

rm -rf "$WORK"
echo "[$(date '+%H:%M:%S')] Done. Outputs in $OUT_DIR"
