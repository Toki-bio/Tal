# rsi_main_20261022

*Rhinolophus sinicus* (rsi), complete SINEderella rerun with a 22-consensus bank: the 21 consensuses of
`rsi_main_20261021b` plus **r10_r8_P18del** (255 bp). Same genome (`GCF_057876585.1`), the same family table
(`rsi_families.tsv`: r# subfamilies and their composites = **Rhin**, the four MEG-* = **MEG**) and the same settings.

Why the new consensus: on the 21b plates of r10_19seqs the "right flank" carried a conserved G-rich, then A-rich stretch.
It was not flank. 535 of the 761 r10 copies were a composite the bank did not have: r10 followed directly by the r8 module
of r10_r8_P18 (P18 positions about 290-423), without P18's middle piece. r10_r8_P18del was built from those copies
(`tools/vendor/sine_consensus.sh`; build directory on the server: `../r10_r8_P18del_build/`).

- Code: SINEderella `main` 4305f5c (2026-10-22).
- Server: therioserver `/home/toki/rhin/rsi_main_20261022/run_20261022_174009` (`../launch.sh`, `--families`, PUBLISH=1,
  THREADS=32; run 17:40-19:28, flank scan stages 1-7 and the report rebuild until 19:48).

## What moved (copies matched between the two runs by coordinate overlap)

| from (21b) | to (this run) | copies |
|---|---|---|
| r10_19seqs | r10_r8_P18del | 359 |
| not assigned (low bitscore) | r10_r8_P18del | 189 |
| r10_r8_P18 | r10_r8_P18del | 131 |
| not assigned (no unanimous vote) | r10_r8_P18del | 36 |
| r6_210seqs | not assigned (no unanimous vote) | 502 |

r10_r8_P18del: 729 assigned copies; 635 SINE, 78 grey, 16 not SINE (flank twins 2.2 %); both sides tested independent
(step 7). r10_19seqs keeps 315 copies (median 149 bp), r10_r8_P18 4,138.

**r6 copies that failed the vote: run-to-run spread, not the new consensus.** 502 copies that were r6_210seqs in 21b
fail the unanimous vote here; their other votes went to r5_r6_P26 (889 copies with an r6/P26 split here, 423 in 21b), not to
P18del. Each consensus is searched on its own, so a new consensus cannot change r6's or P26's scores. Repeating only the
r6 and P26 cycles (10 cycles, 5 times on each run's own copy library; therioserver `/home/toki/tmp/r6vote/`, `vt.sh`):
copies won by r6 in most but not all cycles, 2,800-3,846 on the 21b library and 2,524-3,870 on this run's library. The
spread between repeats on one library is larger than the difference between the runs.

## Verdicts (every copy; `loci_verdict.tsv`, per subfamily `flank_verdict.tsv`)
Rules: SINEderella `docs/VERDICTS.md`.

| | rsi_main_20261021b | this run |
|---|---|---|
| copies in the list | 67,113 | 67,192 |
| assigned | 63,816 | 63,625 |
| SINE loci | 52,100 | 52,017 |
| grey zone | 9,173 | 9,085 |
| not SINE (assigned + not assigned) | 5,840 | 6,090 |

Per family (assigned copies):

| family | subfamilies | assigned | SINE | grey | not SINE |
|---|---|---|---|---|---|
| Rhin | 18 | 62,759 | 51,896 | 8,354 | 2,509 |
| MEG | 4 | 866 | 121 | 731 | 14 |

## Copies that continue past their consensus (`continuation/`)
This run predates the continuation step (merged into SINEderella main as 3f2b036). The check was run on this run's
copies afterwards with that code (therioserver `/home/toki/tmp/cont_oracle2/22`; SINEderella `docs/CONTINUATION.md`):

- the 21b flag (r10 copies continuing into P18 270-415, 58 %) is gone;
- 128 of 315 r10 copies (41 %) continue into P18's middle piece (positions 173-291), entirely outside their loci: these
  loci stop at r10's end. The candidate built from them (`candidates.fa`) is not in the bank;
- 19 of 85 tested unassigned copies with best match r10 continue into P18 232-361.

## r10_19seqs: what the "simple repeats" in the right flank are (2026-10-22)
On the r10 plates the sequence past r10's end looks like simple repeats. It is the element itself. Measured on the 315
assigned r10 copies, padded 1 kb on the 3' side (therioserver `/home/toki/tmp/tsd_r10/`, `tsdscan.py`, `tsd4.py`):
- words of the GC-rich unit that ends the r10 consensus (GCTGC, GAGC, GTGGGC, GGCAGC, GCCGTG) are 18 times denser in the
  first 25 bp past r10's end than in the 5' flank (1.26 against 0.07 per 25 bp), 12 times at 25-49 bp, and fall to
  about twice the flank level by 150 bp: a tract of r10's own terminal units, of variable length;
- target-site duplications (5' flank repeated on the 3' side, 10-20 bp, <= 10 % mismatches, against a shuffled-pair
  null): 68 of 260 copies with a complete r10 match have one within 1 kb; real hits exceed the null from 25 to 225 bp past
  r10's end, only 3 copies have their nearest TSD within 25 bp of it.
So most r10 insertions end 25-225 bp after the r10 consensus: the consensus stops inside the element, and the workflow
takes the rest (the unit tract, then the tail) as flank. The continuation step's "r10 into P18 173-291" is the same tract
(P18's middle piece is built of these units).

## Plates
`alignments/`: top 100, random 100 and subfamily plates (SINE loci, topped up with grey copies marked [grey]).
No copy row has a gap inside a flank. In 91 of 3,960 rows (21b: 83 of 3,770) there is a run of gap columns between the
last element base and the flank (64 on the 3' side, 27 on the 5' side); not yet fixed.

## Contents
| Path | What |
|---|---|
| `report.html` | full report (subfamily table by family, plates, element hierarchy, consensus similarity, tRNA / 5S / 7SL origin, divergence, step1-2 tables) |
| `loci_verdict.tsv`, `flank_verdict.tsv` | per-copy and per-subfamily verdicts |
| `rsi_families.tsv`, `families.summary.tsv` | the family table given to the run, and the family roll-up summary |
| `alignments/` | plates |
| `continuation/` | continuation table, candidates, the copies each was built from |
| `flankscan/` | hierarchy, peaks, candidates, chains, accepted candidates |
| `assignment_stats.tsv`, `summary.by_subfam.tsv` | per-consensus numbers (assigned copies) |
| `flank_twins.tsv`, `array_flag.tsv`, `satellites/` | stage-9 flank twins per family, tandem-array share, satellite screen |
| `length_variants/`, `consensus_audit/` | length-version test; each consensus rebuilt from its own copies |
| `composite_units.tsv`, `consensus_blocks.tsv` | composite consensuses read as units; the bank drawn as elementary blocks |
| `consensuses.clean.fa`, `manifest.txt` | the bank; tool versions and the command |
