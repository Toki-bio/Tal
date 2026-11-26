# Apple (mdo) peel of the AnnoSINE_v2 seed plate, round 1 (2026-11-26)

Plate: `../alignments/mdo_consensuses.aln.fa` (496 AnnoSINE_v2 seed consensuses, 1,311 columns). Unlike the earlier
peels, the rows are de novo seeds, not SubFam chunk consensuses.

## Owner's result

Seven consensuses, verbatim in `peel_round1_owner_consensuses.fa`; each name carries the owner's member count
(3, 4, 7, 17, 25, 34, 354 = 444 of 496). Member lists, call order and comments were not given.

## Reconstructed membership (not the owner's calls)

`seeds_vs_peeled.py` -> `seeds_vs_peeled.tsv`: each seed goes to its best consensus if identity >= 0.6 over >= 0.5 of
the seed. This gives 6, 6, 16, 34, 4, 24, 352 for the groups named 4, 7, 17, 34, 3, 25, 354 (442 seeds; 54 outside).
The counts agree to within 2 per group, which does not prove that the same seeds are in each group.

`peel_compare.js` -> `peel_compare.txt`: per group, identity within the group and against the rest, diagnostic columns
(definition of `SINE_discriminator/RSI_PEEL_LOG.md`: >= 90 % in the group, <= 50 % in the rest) and position on the
plate, then ViewAlign's own methods scored against the reconstructed groups.

| method (ViewAlign code, minimum group 3) | groups | adjusted Rand index (ARI) | unassigned |
|---|---|---|---|
| tree cut, automatic | 2 | 0.000 | 2 |
| tree cut, 7 groups given | 7 | 0.946 | 23 |
| diagnostic-column peel loop | 8 | 0.787 | 12 |

**Do not read the tree cut's 0.946 as reproducing the owner's judgement.** The reference groups were rebuilt by
sequence similarity, which favours a similarity method, and the number of groups (7) was given to it. A fair score
needs the owner's member lists.

## Against the earlier peel logs

| point | earlier peels (ccr, oma, rsi, hla) | this peel |
|---|---|---|
| rows left out | outliers, noise and singletons are deleted | 52 of 496 are in no named group |
| group shape | a block of adjacent rows, taken out as one | the 354, 25, 17 and 34 groups are mostly adjacent runs (see `peel_compare.txt`) |
| smallest group | two single chunks were "not enough to support a subfamily" (hla) | groups of 3 and 4 seeds were kept |
| scores on ccr / oma | one-shot tree cut ARI about 0.42-0.67, peel loop 0.55-0.84 (ViewAlign v237 notes) | order reversed, with the caveat above |

## Relation to the apple Repbase entries (names only; sequences are licensed)

| peeled consensus | closest apple entry | identity | share of the peeled consensus aligned |
|---|---|---|---|
| consensus_354seqs | SINE-1_Mad | 94.7 % | 99 % |
| consensus_17seqs | SINE-4_Mad | 95.8 % | 98 % |
| consensus_4seqs | SINE-5_Mad (opposite strand) | 97.3 % | 83 % |
| consensus_7seqs | SINE-2_Mad | 81.9 % | 100 % |
| consensus_34seqs | SINE-2_Mad | 74.2 % | 93 % |
| consensus_3seqs | no real match (best 49.6 % over half) | | |
| consensus_25seqs | no real match (best 51.1 %) | | |

## SINEderella run

Bank: the seven consensuses plus Repbase SINE-1_Mad and SINE-3_Mad (owner, 2026-11-26); all nine kept by the
canonicalisation step (no merges). Run: `Biotech2024:/staging/tmp/apple_sines/sd_peel9/run_20261126_011507`, SINEderella
b5c8a4a. consensus_354seqs and SINE-1_Mad compete for the same copies.

### Result (finished 2026-11-26, exit 0)

19,335 extracted copies, 15,113 assigned; verdicts: SINE 12,225, grey 1,730, not SINE 5,380.

| query | firm copies | soft | leak % of firm | verdict SINE | consensus audit |
|---|---|---|---|---|---|
| SINE-1_Mad | 8,781 | 1,811 | 34.3 | 7,094 | SHORTER |
| consensus_25seqs | 1,804 | 520 | 0.0 | 1,489 | MATCH |
| consensus_4seqs | 1,781 | 55 | 0.0 | 1,560 | MATCH |
| consensus_7seqs | 1,200 | 128 | 2.4 | 1,022 | DIVERGED |
| consensus_3seqs | 450 | 263 | 0.0 | 112 | LONGER |
| consensus_354seqs | 381 | 951 | 85.3 | 331 | DIVERGED |
| consensus_34seqs | 376 | 38 | 0.5 | 336 | MATCH |
| consensus_17seqs | 324 | 67 | 2.5 | 266 | SHORTER |
| SINE-3_Mad | 16 | 389 | 0.0 | 15 | SKIPPED |

Length-version test: consensus_354seqs against SINE-1_Mad = one end mode (no second length version).
consensus_3seqs: 330 of 400 sampled assigned copies continue 5' into the same consensus (dimer or tandem copies);
323 of its 450 assigned copies have masked (repeat) flanks and could not be tested.

## Run without SINE-1_Mad (2026-11-26, owner: "remove SINE-1_Mad, keep mine")

`SINEderella --exclude SINE-1_Mad` -> `Biotech2024:/staging/tmp/apple_sines/sd_peel9/run_excl_20261126_020256`, then
`publish_run.sh ... mdo`. consensus_354seqs: 9,203 firm copies, leak 0.01 % (before: 381 and 85.3 %; SINE-1_Mad had
8,781). Verdicts: SINE 12,331, grey 1,715, not SINE 4,841 (18,887 copies). Published copy: `../sd/`, made by
`../strip_repbase.py` (Repbase consensus rows and the SINE-3_Mad sequence string removed).

Seen while checking the published copy (SINEderella b5c8a4a):
- the excluded run's all-seeds plate still carried the SINE-1_Mad consensus row;
- the report's satellite table still lists SINE-1_Mad (carried over from the source run);
- the report's similarity matrix and array table cut names at the first "_", so every `consensus_Nseqs` shows as
  "consensus" (`report_hierarchy.py` `_short`).
