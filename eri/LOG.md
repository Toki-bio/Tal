# eri (Erinaceus europaeus) — analysis log

*Scoped to this species only — see the main [WORKFLOW.md](../WORKFLOW.md) for the
shared Tal SINE pipeline this reuses. Genome: GCF_950295315.1 (mEriEur2.1).*

---

## 2026-08-20 — Initial run: hit a real SubFam bug

First `SINEderella` invocation on this genome
(`GCF_950295315.1_mEriEur2.1_genomic.fna`, 2 known consensuses ERI-1/ERI-2 from
Borodulina 2001, run directory `run_20260820_204827` on KIT) completed genome
sanitization, ran Step 1 search & extract successfully (1,548,346 candidate loci
extracted, 251,355 hits for ERI-1 / 1,312,311 hits for ERI-2), then **failed** at
Step 1's SubFam clustering sub-step:

```
/data/V/toki/bin/SINEderella/SubFam: 7: /data/V/toki/bin/SINEderella/SubFam: [[: not found
/data/V/toki/bin/SINEderella/SubFam: 11: /data/V/toki/bin/SINEderella/SubFam: Syntax error: "(" unexpected
```

**Root cause**: `SubFam`'s shebang was `#!/bin/sh`, but the script body used
bash-only syntax (`[[ ]]` test brackets, array assignment). On KIT, `/bin/sh` is
dash, which doesn't support either — 100% reproducible on every invocation.
Confirmed identical bug in the canonical GitHub copy
(`Toki-bio/SINEderella`), not local drift. Also confirmed the same failure had
already silently occurred in an earlier, unrelated run (the dog/*Canis
familiaris* validation run) without being noticed, since that run's Step 2
onward was resumed manually, bypassing the orchestrator's fatal-exit check on
Step 1.

**Fixed**: shebang changed to `#!/bin/bash` (single-line diff). Deployed to the
live shared tool path (`/data/V/toki/bin/SINEderella/SubFam`, PATH-resolved by
every invocation, not a per-run frozen copy) and committed/pushed to
`Toki-bio/SINEderella` (`main`, commit `8feb497`/merge `8c65436`), alongside two
other fixes made the same session (LEAK-detection dead-code fix in
`step2_asSINEment.sh`, python3-version auto-detect fix in `step4_plots.sh`).

A second attempt to restart the full `SINEderella` pipeline was interrupted by
the user mid-genome-sanitization (`run_20260820_210902`) — that run's
`genome.clean.fa` was left at 0 bytes (confirmed directly), so it was abandoned
rather than resumed (SINEderella's `--resume` mode only checks file
existence/non-emptiness, not completeness, so resuming a truncated genome file
would have silently run the full pipeline on corrupt input).

## 2026-08-20 — Reused Step 1's output, ran Step 2 manually

Rather than redo the expensive Step 1 search from scratch, reused the *original*
run's already-complete `genome.clean_step1/extracted.fasta` (1,548,346 sequences,
produced before SubFam killed the run) directly as Step 2's input — same pattern
used for the dog genome validation run. Deployed the fixed (LEAK-tracking)
`step2_asSINEment.sh` into `run_20260820_204827/`, launched:

```bash
bash step2_asSINEment.sh consensuses.clean.fa genome.clean_step1/extracted.fasta 48 step2/step2_output
```

**Result**: 1,548,011/1,548,346 unanimous (10/10 votes), 1,276,402 passing all
criteria (82.4%) — ERI-1: 236,299 sequences, ERI-2: 1,311,712 sequences.

## 2026-08-20 — Ran the now-fixed SubFam directly against the original sample

`genome.clean_step1/subfam_input/input.fasta` (the 30,000-sequence sample SubFam
was originally given before crashing) still existed from the first attempt.
Ran the fixed `SubFam` directly against it (not through the full orchestrator,
since Step 1 as a whole was bypassed):

```bash
export PATH=/data/V/toki/bin/SINEderella:$PATH
cd run_20260820_204827/genome.clean_step1/subfam_input
SubFam input.fasta 50
```

Completed successfully in 4m9s ("Reordering bank and making consensus
sequences" → "Aligning consensus sequences" → "Final alignment" → "Job
completed"), producing `input.clw` — an alignment of ~600 per-chunk consensus
sequences (each chunk MAFFT-aligned independently, then all chunk-consensuses
aligned together in a final cross-chunk step), not a raw alignment of the
original 30,000 sequences and not true CLUSTAL block format despite the `.clw`
extension (plain gapped FASTA). Downloaded locally.

## 2026-08-20/21 — Manual clustering into candidate subfamilies

User manually inspected `input.clw` and clustered the per-chunk consensus
sequences into two candidate subfamily groups:
- **ERI-1-like**: e1-1 (36 seqs, closest match to the literature ERI-1
  consensus), e1-2 (19 seqs), e1-3 (32 seqs), e1-4 (3 seqs)
- **ERI-2-like**: e2-1 (225 seqs), e2-2 (4 seqs), e2-3 (231 seqs), e2-4 (44 seqs)

Result aligned against the original SINEbase/literature ERI-1 and ERI-2
consensus sequences (Borodulina 2001) — see
[`alignments/eri_candidate_consensuses_vs_sinebase.aln.fa`](alignments/eri_candidate_consensuses_vs_sinebase.aln.fa).

**Next step (user, not yet done)**: run `SINEderella` against this new,
expanded candidate consensus set (8 sequences instead of 2) to get a real
subfamily classification, then extract top100/rand100/subfam alignments per
subfamily (see [`scripts/`](scripts/)) and build a full `report.html` matching
the other 6 Tal species.

## 2026-08-21 — Second SINEderella run against the 8-consensus candidate set

Ran the full pipeline (`run_20260820_221537` on KIT) using
`eri_candidates8.fasta` → `consensuses.clean.fa` (8 sequences: e1-1..e1-4,
e2-1..e2-4) as the reference bank. Step 1 (search+extract) and Step 2
(assignment) both completed. All 8 candidates received real genomic hits —
none came back empty:

| Subfamily | Assigned copies |
|---|---|
| e1-1 | 123,749 |
| e1-2 | 146,021 |
| e1-3 *(header has no space: `e1-3consensus_32seqs`, inherited from the source consensus file, not fixed here)* | 83,033 |
| e1-4 | 10,859 |
| e2-1 | 631,692 |
| e2-2 | 11,578 |
| e2-3 | 647,525 |
| e2-4 | 81,431 |

## 2026-08-21 — Extraction script: found and adopted the real Tal precedent

The first version of `scripts/extract_top100_rand100_subfam.sh` (staged
2026-08-20, described above) was written without an exact precedent —
its own header said so explicitly. Two problems surfaced when actually
running it against real data:

1. **Wrong sampling method.** It used a custom seeded Python
   `random.sample()` for `rand100`. The project has its own dedicated tool
   for this, `/data/V/toki/bin/sample N file.fa` (unseeded, `shuf`-based),
   which is the one that should be used for reproducing the same convention
   as the rest of the site.
2. **Missing genomic flanks.** The first version aligned the raw matched
   region only. The real Tal species pages include 50bp-left/70bp-right
   genomic flanks around each copy (for TSD inspection) — confirmed by two
   independent pieces of evidence: the Tal repo's own commit history
   (`"Add flanked alignments (50+70bp)"`, `"...flanked versions (50bp left,
   70bp right)"`) and direct user confirmation.

Searched KIT directly for the actual script that produced the existing
saq/ccr/etc. alignments (git history in the Tal repo itself only ever
contains the *output* `.aln.fa` files, never a generation script — so this
had to be found on the server, not in either git repo). Found
`/data/V/toki/bin/SINEderella/extract_alignments.sh`, confirmed identical
(md5-matched, modulo text encoding) to `/data/W/toki/
extract_alignments_sq.sh` — the actual script invoked for a prior Tal
run (log at `/data/W/toki/extract_alignments_sq.log`). It implements exactly
this: re-extract flanked copies from the genome via `bedtools slop`+
`getfasta` using coordinates parsed from the SINEderella header format
(`>ctg:start-end(strand)|Subfam|Bits`), then a post-alignment cleanup step
that lowercases and re-justifies the flank bases against the consensus body
boundary. Its own header comment says "30bp left", which conflicts with the
50bp confirmed by the commit history/user — used 50, not 30.

`scripts/extract_top100_rand100_subfam.sh` was rewritten to fold in this
exact real logic (functions renamed/kept close to the original for
traceability) instead of the original from-scratch reimplementation. Key
parts, embedded here since they're the actual mechanism behind every
`top100`/`rand100` alignment on this page:

```bash
# Parse FASTA headers to BED. Header format: >ctg:start-end(strand)|Subfam|Bits
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

# Re-extract with FLANK_L=50 / FLANK_R=70 flanks from the genome
extract_flanked(){
  local infile="$1" outfile="$2" tmpd="$3"
  fasta_headers_to_bed "$infile" "$tmpd/flank.bed"
  [[ -s "$tmpd/flank.bed" ]] || { cp "$infile" "$outfile"; return; }
  awk -v OFS='\t' '{print $1,$2}' "${GENOME}.fai" > "$tmpd/genome.sizes"
  bedtools slop -s -l "$FLANK_L" -r "$FLANK_R" -g "$tmpd/genome.sizes" \
    -i "$tmpd/flank.bed" > "$tmpd/flanked.bed" 2>/dev/null || true
  [[ -s "$tmpd/flanked.bed" ]] && bedtools getfasta -s -nameOnly -fi "$GENOME" \
    -bed "$tmpd/flanked.bed" > "$outfile" 2>/dev/null || true
  [[ -s "$outfile" ]] || cp "$infile" "$outfile"
}
```

Full script (including the `postprocess_flanks` alignment-cleanup step, also
carried over from the real precedent): [`scripts/extract_top100_rand100_subfam.sh`](scripts/extract_top100_rand100_subfam.sh).

Re-ran against `run_20260820_221537` (in progress as of this entry) —
`e1-1` completed cleanly with flanked `top100`/`rand100` alignments plus a
`subfam` SubFam re-clustering alignment, confirming the corrected script
works end-to-end on real data. Remaining subfamilies still processing;
outputs will be added to [`alignments/`](alignments/) once complete.
