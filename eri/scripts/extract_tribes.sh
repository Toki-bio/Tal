#!/bin/bash
# extract_tribes.sh — "SINE Tribes" (99% identity copy clusters) extraction,
# matching the real method used for saq (see saq/tribes/tribes_metadata.txt,
# which is the actual, verified spec this script implements — found by
# reading the committed metadata file, not reverse-engineered/guessed):
#
#   method=SINEtribes-compatible vsearch --cluster_fast --id 0.99 --strand plus
#   min_cluster_size=5
#   top_alignment_tribes=10
#   max_alignment_sequences=500
#   sampling=awk/shuf
#   alignment_flanks=50bp_left_70bp_right
#   alignment_flank_method=strand-aware genome.clean.fa extraction; lowercase
#     unaligned flank blocks attached to existing MAFFT body alignment
#
# Unlike top100/rand100/subfam (which cluster/sample from step2's per-
# subfamily *assigned* output), Tribes clusters the RAW Step1 extracted pool
# (genome.clean_step1/extracted.fasta) — i.e. sequence-identity clusters that
# exist independent of subfamily assignment. Each tribe's dominant_subfamily
# is then computed by cross-referencing step2's assigned.fasta, purely as
# descriptive annotation.
#
# Output (matches saq/tribes/ exactly):
#   {species}_tribes_summary.tsv
#   {species}_tribeNN_<n>seqs.aln.fa   (top 10 tribes by size, flanked+aligned)
#   tribes_metadata.txt
#
# Usage:
#   extract_tribes.sh <RUN_ROOT> <SPECIES_CODE> <OUT_DIR>

set -euo pipefail

RUN_ROOT="${1:?usage: $0 <RUN_ROOT> <SPECIES_CODE> <OUT_DIR>}"
SPECIES="${2:?usage: $0 <RUN_ROOT> <SPECIES_CODE> <OUT_DIR>}"
OUT_DIR="${3:?usage: $0 <RUN_ROOT> <SPECIES_CODE> <OUT_DIR>}"

EXTRACTED="$RUN_ROOT/genome.clean_step1/extracted.fasta"
ASSIGNED="$RUN_ROOT/step2/step2_output/assigned.fasta"
GENOME="$RUN_ROOT/genome.clean.fa"
VSEARCH_BIN="/data/V/toki/miniforge3/bin/vsearch"
MIN_CLUSTER_SIZE=5
TOP_TRIBES=10
MAX_ALN_SEQS=500
FLANK_L=50
FLANK_R=70
THREADS="$(nproc)"

[[ -s "$EXTRACTED" ]] || { echo "ERROR: missing/empty $EXTRACTED" >&2; exit 1; }
[[ -s "$ASSIGNED"  ]] || { echo "ERROR: missing/empty $ASSIGNED" >&2; exit 1; }
[[ -s "$GENOME.fai" ]] || { echo "ERROR: missing $GENOME.fai" >&2; exit 1; }
[[ -x "$VSEARCH_BIN" ]] || { echo "ERROR: vsearch not found at $VSEARCH_BIN" >&2; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not in PATH" >&2; exit 1; }

PYTHON3="python3"
command -v "$PYTHON3" >/dev/null 2>&1 || {
  for cand in python3.13 python3.12 python3.11 python3.10 python3.9; do
    command -v "$cand" >/dev/null 2>&1 && { PYTHON3="$cand"; break; }
  done
}

mkdir -p "$OUT_DIR"
WORK="$OUT_DIR/_work"
rm -rf "$WORK"
mkdir -p "$WORK"

# ---- Same flank machinery as extract_top100_rand100_subfam.sh ----
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
    while ((p = index(rest, ":")) > 0) { last_colon += p; rest = substr(rest, p+1) }
    if (last_colon == 0) next
    ctg = substr(loc, 1, last_colon-1)
    coords = substr(loc, last_colon+1)
    split(coords, c, "-")
    if (c[1]+0 >= 0 && c[2]+0 > 0) {
      printf "%s\t%s\t%s\t%s\t0\t%s\n", ctg, c[1], c[2], id, strand
    }
  }' > "$outfile"
}

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
  [[ -s "$outfile" ]] || cp "$infile" "$outfile"
}

echo "[$(date '+%H:%M:%S')] Clustering $EXTRACTED with vsearch --cluster_fast --id 0.99 --strand plus ..."
"$VSEARCH_BIN" --cluster_fast "$EXTRACTED" --id 0.99 --strand plus \
  --threads "$THREADS" --uc "$WORK/clusters.uc" --quiet

echo "[$(date '+%H:%M:%S')] Building seqid->subfamily map from $ASSIGNED ..."
"$PYTHON3" - "$ASSIGNED" "$WORK/seq_to_subfam.tsv" <<'PYEOF'
import sys
assigned, out = sys.argv[1], sys.argv[2]
with open(assigned) as f, open(out, 'w') as fh:
    for line in f:
        if line.startswith('>'):
            parts = line[1:].strip().split('|')
            if len(parts) >= 2:
                fh.write('{}\t{}\n'.format(parts[0], parts[1]))
PYEOF

echo "[$(date '+%H:%M:%S')] Parsing vsearch .uc, ranking clusters, writing top $TOP_TRIBES ..."
"$PYTHON3" - "$WORK/clusters.uc" "$WORK/seq_to_subfam.tsv" "$MIN_CLUSTER_SIZE" "$TOP_TRIBES" "$WORK" <<'PYEOF'
import sys
from collections import defaultdict

uc_file, map_file, min_size, top_n, work = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]

seq_to_sf = {}
with open(map_file) as f:
    for line in f:
        seqid, sf = line.rstrip('\n').split('\t')
        seq_to_sf[seqid] = sf

clusters = defaultdict(list)  # cluster_id -> [seqid, ...]
with open(uc_file) as f:
    for line in f:
        cols = line.rstrip('\n').split('\t')
        rtype = cols[0]
        if rtype not in ('S', 'H'):
            continue
        cluster_id = cols[1]
        seqid = cols[8].split()[0]
        clusters[cluster_id].append(seqid)

ranked = sorted(
    ((cid, members) for cid, members in clusters.items() if len(members) >= min_size),
    key=lambda kv: len(kv[1]), reverse=True
)

with open(work + '/tribes_ranked.tsv', 'w') as fh:
    for rank, (cid, members) in enumerate(ranked[:top_n], start=1):
        counts = defaultdict(int)
        assigned_n = 0
        for seqid in members:
            sf = seq_to_sf.get(seqid)
            if sf:
                counts[sf] += 1
                assigned_n += 1
        seq_count = len(members)
        unassigned_n = seq_count - assigned_n
        if counts:
            dom_sf, dom_n = max(counts.items(), key=lambda kv: kv[1])
        else:
            dom_sf, dom_n = 'NA', 0
        breakdown = ','.join('{}:{}'.format(sf, n) for sf, n in sorted(counts.items(), key=lambda kv: -kv[1]))
        pct_total = round(100.0 * dom_n / seq_count, 2) if seq_count else 0.0
        pct_assigned = round(100.0 * dom_n / assigned_n, 2) if assigned_n else 0.0
        fh.write('\t'.join(str(x) for x in [
            rank, 'tribe{:02d}'.format(rank), cid, seq_count, assigned_n, unassigned_n,
            dom_sf, dom_n, pct_total, pct_assigned, breakdown
        ]) + '\n')
        with open('{}/tribe{:02d}.members.txt'.format(work, rank), 'w') as mf:
            mf.write('\n'.join(members) + '\n')
PYEOF

echo -e "rank\ttribe_id\tcluster_id\tseq_count\tassigned_count\tunassigned_or_missing_count\tdominant_subfamily\tdominant_subfamily_count\tdominant_subfamily_pct_total\tdominant_subfamily_pct_assigned\tsubfamily_breakdown" \
  > "$OUT_DIR/${SPECIES}_tribes_summary.tsv"
cat "$WORK/tribes_ranked.tsv" >> "$OUT_DIR/${SPECIES}_tribes_summary.tsv"

MAFFT_OPTS=(--thread "$THREADS" --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --preservecase --quiet)

while IFS=$'\t' read -r rank tribe_id cluster_id seq_count rest; do
    [[ -n "$rank" ]] || continue
    members_file="$WORK/${tribe_id}.members.txt"
    [[ -s "$members_file" ]] || continue

    "$PYTHON3" - "$EXTRACTED" "$members_file" "$WORK/${tribe_id}.raw.fasta" <<'PYEOF'
import sys
fasta, members_file, out = sys.argv[1], sys.argv[2], sys.argv[3]
wanted = set(l.strip() for l in open(members_file) if l.strip())
keep = False
with open(fasta) as f, open(out, 'w') as fh:
    for line in f:
        if line.startswith('>'):
            seqid = line[1:].strip().split()[0]
            keep = seqid in wanted
        if keep:
            fh.write(line)
PYEOF

    n_extracted="$(grep -c '^>' "$WORK/${tribe_id}.raw.fasta")"
    if (( n_extracted > MAX_ALN_SEQS )); then
        "$PYTHON3" - "$WORK/${tribe_id}.raw.fasta" "$MAX_ALN_SEQS" "$WORK/${tribe_id}.sampled.fasta" <<'PYEOF'
import sys, random
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
sampled = random.sample(recs, n)
with open(out, 'w') as fh:
    for nm, sq in sampled:
        fh.write('>{}\n{}\n'.format(nm, sq))
PYEOF
        aln_input="$WORK/${tribe_id}.sampled.fasta"
        n_used="$MAX_ALN_SEQS"
    else
        aln_input="$WORK/${tribe_id}.raw.fasta"
        n_used="$n_extracted"
    fi

    extract_flanked "$aln_input" "$WORK/${tribe_id}.flanked.fasta" "$WORK"
    mafft "${MAFFT_OPTS[@]}" "$WORK/${tribe_id}.flanked.fasta" \
        > "$OUT_DIR/${SPECIES}_${tribe_id}_${n_used}seqs.aln.fa"
    echo "[$(date '+%H:%M:%S')] $tribe_id: $seq_count total members, $n_used aligned (50L+70R flanked)"
done < "$WORK/tribes_ranked.tsv"

cat > "$OUT_DIR/tribes_metadata.txt" <<META
species=$SPECIES
run_root=$RUN_ROOT
input=$EXTRACTED
method=SINEtribes-compatible vsearch --cluster_fast --id 0.99 --strand plus
min_cluster_size=$MIN_CLUSTER_SIZE
top_alignment_tribes=$TOP_TRIBES
max_alignment_sequences=$MAX_ALN_SEQS
sampling=python random.sample
alignment_flanks=${FLANK_L}bp_left_${FLANK_R}bp_right
alignment_flank_method=strand-aware genome.clean.fa extraction via bedtools slop+getfasta
alignment_flank_source=genome.clean.fa
META

rm -rf "$WORK"
echo "[$(date '+%H:%M:%S')] Done. Outputs in $OUT_DIR"
