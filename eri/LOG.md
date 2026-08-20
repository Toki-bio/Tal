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
