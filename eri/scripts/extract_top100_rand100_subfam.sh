#!/bin/bash
# extract_top100_rand100_subfam.sh
#
# Per-subfamily result extraction for the Tal SINE site, generalized from
# _sinederella/patch/run_subfam_per_sf.sh (which only covers the "subfam"
# variant and is hardcoded to the saq/ccr runs). This script produces all
# three variants used on the existing species pages (saq/ccr/toc/teu/gpy/dmo):
#
#   {species}_{sf}_top100.aln.fa   — 100 highest-bitscore copies (genome-flanked)
#                                     + consensus, mafft-aligned
#   {species}_{sf}_rand100.aln.fa  — 100 randomly-sampled copies (genome-flanked)
#                                     + consensus, mafft-aligned
#   {species}_{sf}_subfam.aln.fa   — SubFam re-clustering of up to N_SAMPLE copies
#                                     (chunk-consensuses + sf-consensus, NOT
#                                     flanked — matches extract_alignments.sh's
#                                     Tier C, which has no flank step either)
#
# DELIBERATELY KEPT SEPARATE from _sinederella/patch/ (shared, used by other
# species) — this copy is scoped to eri/ only, per explicit instruction not to
# let eri-specific work spread into other Tal species pages.
#
# Method (updated 2026-08-21 after finding the real precedent on KIT):
#   - top100: sort by bitscore parsed from assigned.fasta headers
#     (>ctg:start-end(strand)|subfamily|bitscore) and take the top N.
#   - rand100: the project's own /data/V/toki/bin/sample tool (unseeded,
#     shuf-based N-from-FASTA).
#   - Both top100 and rand100 copies are RE-EXTRACTED FROM THE GENOME WITH
#     FLANKS before alignment (50bp left + 70bp right of the matched region,
#     strand-aware) — this step was missing from the first version of this
#     script. fasta_headers_to_bed/extract_flanked/postprocess_flanks below
#     are adapted directly from the verified real precedent script found on
#     KIT: /data/V/toki/bin/SINEderella/extract_alignments.sh, confirmed
#     identical (md5-matched modulo encoding) to /data/W/toki/
#     extract_alignments_sq.sh, the actual script used for a prior Tal run
#     (log: /data/W/toki/extract_alignments_sq.log). That script's own header
#     comment says "30bp left", but the actual historical commits that added
#     these files to the Tal repo ("Add flanked alignments (50+70bp)",
#     "...flanked versions (50bp left, 70bp right)") and explicit user
#     confirmation both say 50bp left — used 50 here, not the script's 30.
#   - postprocess_flanks (also copied from the real script) lowercases flank
#     bases and re-justifies them against the consensus body boundary after
#     alignment, for TSD inspection — same visual convention as other species.
#   - Both top100 and rand100 get the subfamily consensus appended before
#     the final mafft alignment.
#   - mafft parameters (--localpair --maxiterate 1000 --ep 0.123 --nuc
#     --reorder --preservecase) match the real precedent script exactly.
#   - Contig names in headers get "@U@" substituted for "_" upstream by the
#     SINEderella orchestrator itself (SINEderella:214, gsub(/_/,"@U@",hdr) —
#     a reversible sanitization also used by sear/sear_multi, which restore
#     it on their own outputs). This script's outputs never had that restore
#     step, so top100/rand100 headers leaked "@U@" all the way to the
#     published alignments (e.g. NC@U@080165.1 instead of NC_080165.1).
#     Fixed by restoring "_" right before writing each output file.
#
# Usage:
#   extract_top100_rand100_subfam.sh <RUN_ROOT> <SPECIES_CODE> <OUT_DIR>
#
# <RUN_ROOT>      SINEderella run directory (contains consensuses.clean.fa,
#                  genome.clean.fa(+.fai), step2/step2_output/assigned.fasta)
# <SPECIES_CODE>  short code used in output filenames, e.g. "eri"
# <OUT_DIR>       where to write the *.aln.fa files (created if missing)

set -euo pipefail

RUN_ROOT="${1:?usage: $0 <RUN_ROOT> <SPECIES_CODE> <OUT_DIR>}"
SPECIES="${2:?usage: $0 <RUN_ROOT> <SPECIES_CODE> <OUT_DIR>}"
OUT_DIR="${3:?usage: $0 <RUN_ROOT> <SPECIES_CODE> <OUT_DIR>}"

ASSIGNED="$RUN_ROOT/step2/step2_output/assigned.fasta"
CONS_FA="$RUN_ROOT/consensuses.clean.fa"
GENOME="$RUN_ROOT/genome.clean.fa"
SUBFAM_BIN="/data/V/toki/bin/SINEderella/SubFam"
SAMPLE_BIN="/data/V/toki/bin/sample"
BINSIZE=50
N_SAMPLE=10000       # subfam variant: max copies fed to SubFam clustering
MIN_COPIES_SUBFAM=400  # subfam variant: skip subfamilies smaller than this
TOPN=100
RANDN=100
SEED=42
FLANK_L=50
FLANK_R=70
THREADS="$(nproc)"

[[ -s "$ASSIGNED" ]] || { echo "ERROR: missing/empty $ASSIGNED" >&2; exit 1; }
[[ -s "$CONS_FA"  ]] || { echo "ERROR: missing/empty $CONS_FA" >&2; exit 1; }
[[ -s "$GENOME.fai" ]] || { echo "ERROR: missing $GENOME.fai (samtools faidx it first)" >&2; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not in PATH" >&2; exit 1; }

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

# ---- Parse FASTA headers to BED (adapted from extract_alignments.sh) ----
# Header format: >ctg:start-end(strand)|Subfam|Bits
fasta_headers_to_bed(){
  local infile="$1" outfile="$2"
  grep "^>" "$infile" | sed 's/^>//' | awk '{
    id = $1
    n = split(id, parts, "|")
    loc = parts[1]

    strand = "+"
    if (sub(/\([+-]\)$/, "", loc)) {
      i = index(parts[1], "(")
      strand = substr(parts[1], i+1, 1)
    }

    last_colon = 0; rest = loc
    while ((p = index(rest, ":")) > 0) {
      last_colon += p; rest = substr(rest, p+1)
    }
    if (last_colon == 0) next

    ctg = substr(loc, 1, last_colon-1)
    coords = substr(loc, last_colon+1)
    split(coords, c, "-")
    if (c[1]+0 >= 0 && c[2]+0 > 0) {
      printf "%s\t%s\t%s\t%s\t0\t%s\n", ctg, c[1], c[2], id, strand
    }
  }' > "$outfile"
}

# ---- Re-extract sequences with FLANK_L/FLANK_R flanks from the genome ----
# (adapted from extract_alignments.sh's extract_flanked, flank sizes changed
# from that script's 30/70 to 50/70 per the real Tal-repo commit history)
extract_flanked(){
  local infile="$1" outfile="$2" tmpd="$3"

  fasta_headers_to_bed "$infile" "$tmpd/flank.bed"
  if [[ ! -s "$tmpd/flank.bed" ]]; then
    cp "$infile" "$outfile"
    return
  fi

  awk -v OFS='\t' '{print $1,$2}' "${GENOME}.fai" > "$tmpd/genome.sizes"

  bedtools slop -s -l "$FLANK_L" -r "$FLANK_R" -g "$tmpd/genome.sizes" \
    -i "$tmpd/flank.bed" > "$tmpd/flanked.bed" 2>/dev/null || true

  if [[ -s "$tmpd/flanked.bed" ]]; then
    bedtools getfasta -s -nameOnly -fi "$GENOME" \
      -bed "$tmpd/flanked.bed" > "$outfile" 2>/dev/null || true
  fi

  if [[ ! -s "$outfile" ]]; then
    cp "$infile" "$outfile"
  fi
}

# ---- Post-alignment flank cleanup (verbatim from extract_alignments.sh) ----
# Lowercases flank bases and right/left-justifies them against the consensus
# body boundary, collapsing internal gaps, without changing column count.
postprocess_flanks(){
  local aln_file="$1" cons_name="$2"
  [[ -s "$aln_file" ]] || return 0

  local tmp="${aln_file}.ppf"
  awk -v cons_name="$cons_name" '
    /^>/ {
      if (NR>1) seqs[n] = seq
      n++; hdr[n] = $0; seq = ""
      h = $0; sub(/^>/, "", h); sub(/[[:space:]].*/, "", h)
      if (h == cons_name) cons_idx = n
      next
    }
    { seq = seq $0 }
    END {
      seqs[n] = seq
      if (cons_idx == 0) {
        for (i=1; i<=n; i++) {
          print hdr[i]
          s = seqs[i]; L = length(s)
          for (j=1; j<=L; j+=80) print substr(s, j, 80)
        }
        exit
      }

      cs = seqs[cons_idx]
      alen = length(cs)

      lbound = 0; rbound = 0
      for (j=1; j<=alen; j++) {
        c = substr(cs, j, 1)
        if (c != "-" && c != ".") { if (!lbound) lbound = j; rbound = j }
      }
      if (lbound == 0) lbound = 1
      if (rbound == 0) rbound = alen

      for (i=1; i<=n; i++) {
        print hdr[i]
        s = seqs[i]

        if (i == cons_idx || alen == 0) {
          for (j=1; j<=alen; j+=80) print substr(s, j, 80)
          continue
        }

        lf_bases = ""
        lf_len = lbound - 1
        for (j=1; j<lbound; j++) {
          c = substr(s, j, 1)
          if (c != "-" && c != ".") lf_bases = lf_bases tolower(c)
        }
        lf_pad = lf_len - length(lf_bases)
        lf_out = ""
        for (j=1; j<=lf_pad; j++) lf_out = lf_out "-"
        lf_out = lf_out lf_bases

        body = substr(s, lbound, rbound - lbound + 1)

        rf_bases = ""
        rf_len = alen - rbound
        for (j=rbound+1; j<=alen; j++) {
          c = substr(s, j, 1)
          if (c != "-" && c != ".") rf_bases = rf_bases tolower(c)
        }
        rf_pad = rf_len - length(rf_bases)
        rf_out = rf_bases
        for (j=1; j<=rf_pad; j++) rf_out = rf_out "-"

        full = lf_out body rf_out
        for (j=1; j<=length(full); j+=80) print substr(full, j, 80)
      }
    }
  ' "$aln_file" > "$tmp" && mv "$tmp" "$aln_file"
}

echo "[$(date '+%H:%M:%S')] Splitting $ASSIGNED by subfamily, ranking by bitscore..."

# ---- Split assigned.fasta into per-subfamily FASTA, sorted by bitscore desc ----
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
    cons_name="${sf}_CONSENSUS"
    n_total="$(grep -c '^>' "$sorted_fa")"

    if [[ ! -f "$cons" ]]; then
        echo "[$(date '+%H:%M:%S')] WARN  $sf: no consensus fasta found, skipping"
        continue
    fi

    echo "[$(date '+%H:%M:%S')] $sf: $n_total total assigned copies"

    # ---- top100 (flanked) ----
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
    extract_flanked "$WORK/${sf}.top${TOPN}.fasta" "$WORK/${sf}.top${TOPN}.flanked.fasta" "$WORK"
    cat "$WORK/${sf}.top${TOPN}.flanked.fasta" "$cons" > "$WORK/${sf}.top${TOPN}.combined.fasta"
    mafft "${MAFFT_OPTS[@]}" "$WORK/${sf}.top${TOPN}.combined.fasta" \
        > "$OUT_DIR/${SPECIES}_${sf}_top100.aln.fa"
    postprocess_flanks "$OUT_DIR/${SPECIES}_${sf}_top100.aln.fa" "$cons_name"
    sed -i 's/@U@/_/g' "$OUT_DIR/${SPECIES}_${sf}_top100.aln.fa"
    echo "[$(date '+%H:%M:%S')]   OK top100 ($head_n copies, ${FLANK_L}L+${FLANK_R}R flanked, + consensus)"

    # ---- rand100 (flanked; project's own sample tool: unseeded shuf-based) ----
    rand_n=$(( n_total < RANDN ? n_total : RANDN ))
    cp "$sorted_fa" "$WORK/${sf}.forsample.fasta"
    "$SAMPLE_BIN" "$rand_n" "$WORK/${sf}.forsample.fasta"
    mv "$WORK/${sf}.forsample.fasta.${rand_n}" "$WORK/${sf}.rand${RANDN}.fasta"
    extract_flanked "$WORK/${sf}.rand${RANDN}.fasta" "$WORK/${sf}.rand${RANDN}.flanked.fasta" "$WORK"
    cat "$WORK/${sf}.rand${RANDN}.flanked.fasta" "$cons" > "$WORK/${sf}.rand${RANDN}.combined.fasta"
    mafft "${MAFFT_OPTS[@]}" "$WORK/${sf}.rand${RANDN}.combined.fasta" \
        > "$OUT_DIR/${SPECIES}_${sf}_rand100.aln.fa"
    postprocess_flanks "$OUT_DIR/${SPECIES}_${sf}_rand100.aln.fa" "$cons_name"
    sed -i 's/@U@/_/g' "$OUT_DIR/${SPECIES}_${sf}_rand100.aln.fa"
    echo "[$(date '+%H:%M:%S')]   OK rand100 ($rand_n copies, ${FLANK_L}L+${FLANK_R}R flanked, + consensus)"

    # ---- subfam (SubFam re-clustering, same method as run_subfam_per_sf.sh —
    #      NOT flanked, matches extract_alignments.sh's Tier C) ----
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

# Intermediates ($WORK) are NOT auto-deleted. Clean up manually if truly
# done with them.
echo "[$(date '+%H:%M:%S')] Done. Outputs in $OUT_DIR (intermediates kept in $WORK)"
