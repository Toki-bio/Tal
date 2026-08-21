#!/usr/bin/env bash
# run_trf_genome.sh — Genome-wide Tandem Repeats Finder, parallelized by
# contig (splits into N size-balanced chunk FASTAs, runs trf per chunk via
# GNU parallel, converts to BED, merges). Written to identify tandem-array
# intervals for excluding false-positive "SINE" hits that are actually
# matches to a conserved sub-motif recurring inside a tandem repeat unit
# (see eri/LOG.md 2026-08-21: several Tribes showed near-perfect periodic
# spacing of ~3,425bp and clean multiples of it, not the random spacing
# expected of genuinely dispersed SINE insertions).
#
# MaxPeriod=2000 (standard TRF default), NOT the ~3,425-7,850bp range of
# the macro-periodicity itself — the installed trf 4.09 build segfaults
# above ~2000 (tested: 2000 works, 2200+ crashes consistently, confirmed
# on this exact binary at /data/V/toki/bin/trf). Verified this is still
# sufficient: direct test on a known periodic region (tribe10, ~3,425bp
# SINE-hit spacing) found real tandem content at MaxPeriod=2000 — small
# microsatellite-like periods (1, 3, 4, 10bp) recurring in an alternating
# ~1,160bp/~2,260bp pattern that sums to the same ~3,425bp macro-period.
# The true tandem repeat unit is tiny; the larger inter-hit spacing is
# just the SINE-query-motif's recurrence interval within the array, not
# itself the tandem period TRF needs to detect.
#
# Usage: run_trf_genome.sh <genome.fa> <outdir> [threads]

set -euo pipefail
GENOME="${1:?usage: $0 genome.fa outdir [threads]}"
OUTDIR="${2:?usage: $0 genome.fa outdir [threads]}"
THREADS="${3:-16}"

# Must be absolute: each parallel job does `cd "$OUTDIR/dat"` before using
# a path built from $OUTDIR — a relative OUTDIR breaks after the cd.
GENOME="$(readlink -f "$GENOME")"
OUTDIR="$(mkdir -p "$OUTDIR" && readlink -f "$OUTDIR")"

TRF_BIN=/data/V/toki/bin/trf
TRF2BED=/data/V/toki/bin/trf2bed
MATCH=2 MISMATCH=7 DELTA=7 PM=80 PI=10 MINSCORE=50 MAXPERIOD=2000

mkdir -p "$OUTDIR/chunks" "$OUTDIR/dat" "$OUTDIR/bed"
FAI="${GENOME}.fai"
[[ -s "$FAI" ]] || { echo "ERROR: missing $FAI" >&2; exit 1; }

log(){ printf '[%s] %s\n' "$(date '+%H:%M:%S')" "$*"; }

log "Splitting $GENOME into $THREADS size-balanced chunks by contig..."
python3 - "$FAI" "$THREADS" "$OUTDIR/chunks/assignment.tsv" <<'PYEOF'
import sys
fai, n_chunks, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
contigs = []
with open(fai) as f:
    for line in f:
        cols = line.rstrip('\n').split('\t')
        contigs.append((cols[0], int(cols[1])))
contigs.sort(key=lambda t: -t[1])
bins = [0]*n_chunks
assign = []
for name, length in contigs:
    i = min(range(n_chunks), key=lambda i: bins[i])
    bins[i] += length
    assign.append((name, i))
with open(out, 'w') as f:
    for name, i in assign:
        f.write(f'{name}\t{i}\n')
print('Bin sizes (bp):', bins)
PYEOF

python3 - "$GENOME" "$OUTDIR/chunks/assignment.tsv" "$OUTDIR/chunks" <<'PYEOF'
import sys
genome, assign_file, chunks_dir = sys.argv[1], sys.argv[2], sys.argv[3]
assign = {}
with open(assign_file) as f:
    for line in f:
        name, i = line.rstrip('\n').split('\t')
        assign[name] = int(i)
handles = {}
def get_handle(i):
    if i not in handles:
        handles[i] = open(f'{chunks_dir}/chunk_{i:03d}.fa', 'w')
    return handles[i]
cur = None
with open(genome) as f:
    for line in f:
        if line.startswith('>'):
            name = line[1:].split()[0]
            cur = get_handle(assign[name])
        cur.write(line)
for h in handles.values():
    h.close()
print(f'Wrote {len(handles)} chunk files')
PYEOF

log "Running trf on each chunk (parallel -j $THREADS)..."
# trf exits 1 on SUCCESS (confirmed empirically on this build), so with
# set -e a clean run of every chunk would still abort the script via
# parallel's own nonzero exit — || true is required here, not optional.
ls "$OUTDIR/chunks"/chunk_*.fa | \
  parallel -j "$THREADS" --eta \
    "cd '$OUTDIR/dat' && '$TRF_BIN' {} $MATCH $MISMATCH $DELTA $PM $PI $MINSCORE $MAXPERIOD -d -h > {/}.trf.log 2>&1" \
  || true

log "Converting .dat outputs to BED and merging..."
# trf2bed writes its own ${input%.*}.bed file (named after the .dat, NOT
# printed to stdout despite looking like it should be) — cat the real
# per-chunk .bed files it produces, don't try to capture its stdout.
: > "$OUTDIR/tandem_repeats.bed"
for dat in "$OUTDIR/dat"/chunk_*.fa.*.dat; do
  [[ -s "$dat" ]] || continue
  "$TRF2BED" "$dat" >/dev/null 2>&1 || true
  bed="${dat%.*}.bed"
  [[ -s "$bed" ]] && cat "$bed" >> "$OUTDIR/tandem_repeats.bed"
done
sort -k1,1 -k2,2n "$OUTDIR/tandem_repeats.bed" -o "$OUTDIR/tandem_repeats.sorted.bed"
bedtools merge -i "$OUTDIR/tandem_repeats.sorted.bed" > "$OUTDIR/tandem_repeats.merged.bed"

n_intervals=$(wc -l < "$OUTDIR/tandem_repeats.merged.bed")
total_bp=$(awk '{s+=$3-$2} END{print s+0}' "$OUTDIR/tandem_repeats.merged.bed")
log "Done. $n_intervals merged tandem-repeat intervals, ${total_bp} bp total."
log "Output: $OUTDIR/tandem_repeats.merged.bed"
