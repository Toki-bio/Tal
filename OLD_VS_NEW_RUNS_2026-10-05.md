# Tal page cases: the published runs versus the current SINEderella (2026-10-05)

Written 2026-10-05 for the question "old versus new runs of the Tal page cases, especially bats and scorpions".
Numbers come from the run directories named below and from `SINEderella/docs/RELEASE_CHECK_2026-10-05.md`;
nothing here is estimated. Reruns launched today are listed in section 4 and their results are added to
section 5 as they finish.

## 1. What the pages were built with

| Page group | Species | Report generated | Code state of the run |
|---|---|---|---|
| Talpidae | saq, gpy, toc, dmo, teu | 2026-04-26/27 | April 2026 orchestrator (KIT runs) |
| Talpidae | ccr | 2026-05-18 | May 2026 |
| eri, sco (scorpion *Belisarius xambeui*) | eri, sco | 2026-08-21 | August 2026; sco = 29 hand-clustered candidates, DRAGEN run `run_20260821_025634` |
| zeb, tim, timb, hum | zeb, tim, timb, hum | 2026-08-23/27 | August 2026 (DRAGEN benchmark runs) |
| Chiroptera | 24 bats + rsi | full runs 2026-09-26/27, `--add MEG` 2026-09-27, republished 2026-09-28 18:32 | commits `dd552e5`..`517e041` (09-26/27); publish with `d730836` (installed step6, rand100 seed 42) |
| rsi rebuilds | rsi_v2..v7, rsi_fresh, rsi_fresh2, rsi_final, rsi_sat | 2026-09-29 to 10-04 | up to the satellite screen (rsi_sat, 10-03) |

The bat inputs (therioserver): genomes `~/chiro/genomes/<sp>.fna` (21 species), `~/rhin/genomes/{rsi,rda,rre}.fna`,
`~/lly/genomes/lly.fna`; bank `~/chiro/chiro_bank.fa` (Rhin-1, VES), then `--add ~/chiro/MEG.resolved.fa`
(MEG-RL, MEG-RS, MEG-T2, MEG-TR). rda and rre were run with Rhin-1 + the four MEG only (no VES); lly with all six.
Published run per species: `~/chiro/<sp>/run_add_20260927_*` (rsi: `~/rhin/rsi*`).

## 2. What changed in SINEderella since those runs

All commits are in `Toki-bio/SINEderella`; the current HEAD is `37e144f` (2026-10-05).

| Date | Change | Effect on a rerun |
|---|---|---|
| 09-28 | rand100 drawn with a fixed seed (42); step8a plates parallel; `[soft]`/`[array]` marks; original consensus kept as row 2, changes as proposals | plates only (already in the 09-28 republish) |
| 10-01 | flank scan stages 5-9 (separate tool, not run by the orchestrator) | none on the pages |
| 10-02 | `mafft --threadit 0` on every iterative call | plates reproducible; old plates are not bit-reproducible |
| 10-02 | bank canonicalization fixed: RC pairs merge at >= 80 %, same-orientation only at >= 98 %, length ratio >= 0.9 | the rsi bank lost false merges (r9 into r7, r4 into r2, MEG-RS into MEG-RL); bats: no merge involved the 6-consensus bank |
| 10-02 | length-version test (`tools/length_variants_run.py`) and consensus audit (`tools/consensus_audit.py`) after every assignment | two new tables in the report |
| 10-02 | tandem arrays by regular spacing (`array_order.py`), family flag `array_flag.py` | top-100 takes independent copies first; report "Tandem array" column |
| 10-03/04 | satellite stage in step 1: kind A (TRF unit vs consensus), kind B (regular spacing + unit identity); verified arrays removed from the hit set before extraction and assignment | copy numbers of array-forming families drop; new "Satellites" table |
| 10-05 | audit fixes D1-D27: `satellite_kindB_verify.py` copied into run dirs (the stage failed silently in every orchestrated run from 10-04 until this fix), `--exclude-b verified` default, report columns, `--resume` no longer deletes publish output, array flag with a chance null (dense dispersed families such as VES form regular runs by chance) | correctness of the stage and of the report |

So a current run differs from a 09-27 run in (a) the satellite screen before assignment, (b) three new tables
(satellites, length versions, consensus audit) and the array flag, (c) reproducible plates, and (d) nothing
else in search or assignment: the `-z 11` voting is unchanged, and counts move only within its noise (about 1 %).

## 3. Measured on three bats (release check, 2026-10-04/05)

Fresh clone, same genome, same bank, same `--add MEG.resolved.fa`; runs under `~/tmp/release_check/tal/`.

### rsi (*Rhinolophus sinicus*), 21-consensus bank; reference `~/rhin/rsi_fresh2/run_20261002_112935`

| family | new | reference | note |
|---|---|---|---|
| MEG-RS | 129 | 1,715 | 21 regularly spaced arrays verified (SAT_B) and removed before assignment; the plates' 1,592 array copies are gone |
| r6 | +4.6 % | | shares copies with the removed MEG-RS/MEG-RL loci |
| r5_r5_P48 | 498 | 470 | same |
| all other 18 | within +/- 1 % | | `-z 11` noise |

Length-version verdicts (7 pairs) and consensus-audit verdicts (21 families) identical to the reference.
Satellite table reproduces rsi_sat: MEG-RS 21/21 arrays verified, MEG-TR 1/1, r9 15 small arrays, r4 2.

### tbr (*Tadarida brasiliensis*) and rle (*Rousettus leschenaultii*); references `~/chiro/{tbr,rle}/run_add_20260927_143901`

| family | tbr new / old | rle new / old |
|---|---|---|
| VES | 621,052 / 621,128 | 14 / 14 |
| Rhin-1 | 416 / 437 | 9 / 9 |
| MEG-RL | 10 / 9 | 9,688 / 9,710 |
| MEG-RS | **66 / 249** | **5,949 / 6,718** |
| MEG-T2 | 90 / 92 | 2,153 / 2,159 |
| MEG-TR | 24 / 24 | 2,721 / 2,703 |
| satellites | MEG-RS SAT_B, 2 of 2 arrays verified (81- and 109-unit ~910 bp arrays); MEG-T2 15 kind-A loci; 2,213 hits excluded | MEG-RS 26 of 35 runs verified, MEG-TR 28 of 28; 1,850 hits excluded |
| array flag with chance null | nothing (VES 56.2 % in regular runs vs 47.4 % by chance) | nothing |
| length versions | MEG-RS/MEG-RL NOT_TESTED (too few copies) | MEG-RS/MEG-RL TWO_VERSIONS (Gogolevsky 2009 control reproduced) |
| consensus audit | RS, T2 SHORTER; TR UNSTABLE (24 copies); Rhin-1, VES DIVERGED | RL, RS, TR MATCH; T2 DIVERGED |
| run time (16 threads) | 54 min full + 9 min add | 4 min full + 8 min add |

Reading: the only material change on a bat page is the MEG-RS count where MEG-RS forms arrays (tbr: a quarter
of the old count remains; rle: 11 % fewer). The pages' "Rhin-1 / VES / MEG" numbers for everything else stand.
On the Chiroptera page the MEG-RS rows of species with arrays (tbr, nle, tni, vmu, cse, mme, fho are the
candidates by count) will need the new numbers and the Satellites verdict; which of them actually carry arrays
is what the reruns in section 4 measure.

## 4. Reruns launched 2026-10-05 (therioserver)

Driver `~/tmp/rc_bat_rerun.sh` (detached, log `~/tmp/release_check/bats/driver.log`, marker `bats/DONE`), fresh
clone `~/tmp/release_check/SD` at `13093c9`, 3 species at a time x 16 threads. (A first launch at 07:58 on `37e144f`
was stopped after 10 min and restarted at 08:09 when `13093c9` landed; its partial directories are parked in
`bats/_aborted_37e144f/`.) Per species: full run with
`chiro_bank.fa`, then `--add MEG.resolved.fa`, no publish; `bats/<sp>/summary.txt` holds the new-vs-old table,
satellites, array flag, length versions, audit and warnings.

Species (19): cse cth fho hla mau mev mgi mme msc mtu nle nth ntu rmi rna tni tpe ttr vmu. Started 07:58
server time; the old tbr full run took 20 min on this machine, the new tbr 54 min with the satellite stage, so
expect 4 to 7 hours for the set.

Not in this batch: lly, rda, rre (genomes and old runs under `~/lly/` and `~/rhin/`, different bank
composition; to be run next with their own banks), and rsi/tbr/rle (done, section 3).

### 4b. The 7 kb-period array that the screen missed (other session, 2026-10-05, SINEderella `13093c9`)

On the rsi_sat MEG-RS top-100 plate (`rsi_sat/alignments/rsi_MEG-RS_top100.aln.fa`) 50 of the 100 rows are
units of one array at NC_142507.1:29.82-30.65 Mb with a period of 7,022 bp: flanks 98 % identical between rows,
core 99 %; the other 50 rows are dispersed copies with unrelated flanks (46-49 %) and 68 % core identity. The
regular-spacing rule capped gaps at 6 kb, so the array was never a kind-B candidate. `13093c9` adds a second
tier (runs of >= 10 copies with gaps up to 30 kb replace the narrow runs they contain), used by the screen and
by its chance null; the unit-identity check still decides. Not yet rerun on rsi_sat when this was written; the
reruns in section 4 use it. The same session noted that the array rows carry ~88 bp left flanks against 50 bp
for the other rows and did not trace why (a plate-construction property, open).

Two call sites were not switched by `13093c9`: the plate ordering (`array_order.py` main) and the family flag
(`array_flag.py`) still use the 6 kb rule, so a 7 kb array that is not excluded as "verified" is neither
marked `[array]` by spacing nor counted in the flag. Fixed in SINEderella `26d0cca` (both call `regular_runs_wide`;
18 array tests pass). The three species that started before the pull (cse, cth, fho) get `array_flag.py` re-run
on their finished runs.

## 5. Focus case: MEG-RS in *R. sinicus* (his decision 2026-10-05: one case fixed completely before the next)

The 19-species reruns of section 4 were stopped after 10 minutes (parked in `bats/_stopped_focus_2026-10-05/`). Work
moved to one case. State after the first measurements (yesterday's run `~/tmp/release_check/tal/rsi/run_20261004_233045`,
129 MEG-RS copies left after the 21 verified arrays were removed; stage 9 run on them in `tal/rsi/megrs_fs9/`):

| what | number | meaning |
|---|---|---|
| copies inside a verified array span, not removed | 6 | the exclusion removes only the run members; a hit inside the span that is off the regular spacing stays |
| copies 2-35 kb from a verified array | 27 | array fringes: 50 of them are the "7 kb array" (NC_142507.1:30.57-30.65 Mb, right after the 1 936 bp-unit array ending at 30.56 Mb); five lie 30.6-31.2 kb beyond an array end on five different contigs (a recurring ~31 kb higher-order unit), just past the new 30 kb cap |
| copies > 100 kb from any array or on contigs without one | 82 | the dispersed candidates |
| stage 9 twins (tier 1) | 4 | one group on four contigs, flanks 99-100 % over 93-100 bp, unmasked sequence |
| stage 9 "unique" with the whole 100 bp flank soft-masked on both sides | 11 of the 16 plate rows that share flanks at 96-100 % | fs9 turns soft-masked bases into N and then has nothing to compare, and reports "unique" instead of "untestable"; 79 % of the rsi assembly is soft-masked |
| stage 9 k-mer repeat filter | off in every run so far | `fs9_twins.sh` split the jellyfish dump on tabs, the dump is space-separated; fixed in `5cc42a0` |

Plates (numeric check, `plate_flank_check.py`): rsi_sat top-100: 50 array rows, 16 of the 50 others with a flank partner at
>= 0.96; the release-check top-100: 50 array rows (the 7 kb array), 10 of the 50 others; rand-100: 67 array rows, 9 of 33 others.

Fixes made from this (SINEderella commits of 2026-10-05, each toy-tested on the server before use):

| commit | change | toy check |
|---|---|---|
| `5cc42a0` | `fs9_twins.sh` reads the jellyfish dump split on whitespace (the tab split had left the repeat mask empty in every run) | 4-line dump parsed |
| `7f32adb` | the satellite stage writes `satellites/exclude_regions.bed` and step 1 removes every merged locus inside a kind-A locus or a verified array, whatever consensus found it (the per-consensus filter let array units re-enter through other consensuses' hits); `fs9_twins.sh` keeps soft-masked bases testable by default (SOFTMASK=1 restores) and reports copies with too few testable flank bases as `masked`, not `unique` | planted B1 array gone from the extracted loci, dispersed copies kept; fs9 toy suite passes |
| `c9df4c2` | `flankscan/fs9_run.sh`: stage 9 over every family after assignment, called by the orchestrator (`results/flank_twins.tsv`, per-family `copy_status.tsv`; jellyfish 20-mer mask for genomes > 300 Mb, jellyfish added to the sinederella env); step8a marks twin copies `[twin]` and takes one copy per twin group first | planted segmental duplication (3 copies x 2) = exactly 6 twins, 0 false among 217; marks appear on the toy plates |
| `89607ac` | report table "Are the flanks of the copies unique?" | rendered |

Stage 9 on the 129 old-run copies with the soft-mask lifted: 109 twins (84.5 %), 2 masked, 16 unique. Of the 109, 101 have a
same-contig partner within 50 kb (the 7 kb array and the fringe copies) and 8 are cross-contig: the 4-copy duplicate group and
copies sharing 5' sequence at 85-90 % with array units (old satellite-unit remnants). So after the arrays are gone, the twin
check is what remains to separate duplicates from insertions, and it now runs on every family.

### 5b. MEG-RS mini genome (11.4 Mb of real rsi sequence: +-50 kb around every MEG-RS hit; `therioserver ~/tmp/megrs_mini/`)

Built so that each fix is checked in minutes instead of 85 min on the whole genome (a full run with publish takes 9 + 2.5 min).
What it found:

1. **The 30 kb array tier (13093c9) let a real array back in.** It joined a 2 168 bp-unit MEG-RS array (89 % alone) with the
   dispersed copies 2-35 kb beyond it into one 144-hit run at 82 %, verdict COPIES, so the array was not excluded and filled the
   plates. Fixed in `ed86b68` (narrow runs verified on their own too): the array is now excluded (134 hits, 89 %), MEG-RS keeps
   130 of 1 756 hits after the screen, and 143 of the 236 split-vote loci of the first mini run were in that array.
2. **The same copies get a different vote in a small run.** Of the 129 copies that are firm MEG-RS in the whole-genome run, on the
   mini 89 are MEG-RL with 5 of 10 votes, 6 MEG-RS with 6-7 of 10, 34 not extracted (outside the windows' merged loci); none firm.
   Step 2 searches the consensuses against the copies in library parts of 20 000; with `-z 11` the significance statistics are
   estimated from shuffles of the library, i.e. from 20 000 copies on the whole genome and 351 on the mini. MEG-RS and MEG-RL are
   97 % identical, so with the noisier estimate their order flips between cycles. The strength of the unanimity filter therefore
   depends on how many copies share a library part: the last part of every run, small genomes and the `--add` re-vote (only the
   affected loci) all vote with a weaker or different filter. Not changed: this is the core filter; decision needed.
3. MEG-RS plates on the mini are mostly soft copies (5 firm), so the twin check (firm copies, >= 20) did not cover them: 0 `[twin]`
   marks while 15 of 55 non-array rows share flanks.

The focused rsi run was restarted a second time at 08:53 on `7f32adb` (the 08:29 run lacked the region exclusion); it is read
with the same checks when it finishes, and its MEG-RS plates with `[twin]` marks are the material for the call.

## 6. Scorpions

Three different things exist, none run with the current code:

1. **sco page = *Belisarius xambeui* (GCA_982267015.1)**, 29 hand-clustered candidate consensuses from his SubFam
   run on 22,335 de novo loci; classification run `run_20260821_025634` on DRAGEN (`/staging/tmp/scorpion_sines/`;
   `sco/LOG.md` records that the exact step2 output path was never re-located). 169,262 assigned copies of 251,164
   extracted loci; his manual verdict 2026-08-21: 13 of 29 show any SINE-like signature. August 2026 code: no
   canonicalization, no satellite screen, no length-version or audit tables, plates not reproducible.
2. **oma = *Olivierus martensii* (GCA_000484575.1)**, the manuscript's third demonstration: 23 deduplicated
   AnnoSINE_v2 seeds, 133,843 raw hits, 77,695 merged, 30,000 sampled, 598 chunks
   (`SINE_discriminator/OMA_PEEL_LOG.md`); run `/staging/tmp/scorpions/oma/run_oma` on DRAGEN, repaired 2026-09-08
   (RC-pair merge of the +/- seed duplicates, consensus rebuild, `step4 -3`; `SINEderella/docs/OMA_CONSENSUS_REPAIR.md`
   and `docs/launch_oma_*.sh`). September code.
3. **Nine-genome de novo chain, 2026-10-01** (`/staging/tmp/scorpions_denovo/<code>/`, tables and chunk alignments
   in `sco/denovo_dragen/`): scan + AnnoSINE seeds -> k-mer thinning -> SubFam chunk consensi. Stops before
   assignment; nothing to compare with a SINEderella run.

Why a rerun matters here more than for the bats: scorpion genomes are repeat-rich and the sco candidates were
accepted at 67 % with 16 of 29 judged non-SINE by eye; the satellite screen and the consensus audit are exactly
the tables that would separate array units and unsupported consensuses from SINE families, and the 10-02
canonicalization rules change which of the +/- seed pairs merge (the oma repair used RC >= 80 %, which is still
the rule; same-orientation merges now need >= 98 %).

What blocks it today: DRAGEN did not answer (`plink copilot@100.104.25.22`: "Network error: Connection timed
out", 08:03 local). The genomes, the old runs and a SINEderella checkout (`/staging/tmp/SINEderella`, state
unknown) are there; the same genomes are also on KIT (`/data/V/toki/Genomes/Scorpions/<code>/genome.fna`), where
no SINEderella dependencies have been checked. The current code with every dependency verified is on
therioserver, which has no scorpion genome.

Options, in order of least work: (a) retry DRAGEN later, `git pull` its checkout, verify `ssearch36 seqkit cons
trf samtools gawk`, rerun oma with its 23-seed bank and bxa with the 29-consensus bank, compare with the two old
runs; (b) copy the oma genome (0.9 GB) and the two banks from KIT to therioserver through the KIT tunnel and run
there (about 1 h for oma on 16 threads at this genome size; bxa is 4 GB and would take several hours). Not
launched: both need a server I have not reached today or a 1-5 GB cross-server copy; your call which.

## 7. The other Tal cases

Talpids (saq, gpy, toc, dmo, teu, ccr), eri, zeb, tim/timb and hum were run in April to August 2026 on code
older than every change in section 2. Their runs are on KIT (saq, eri; `SINE_discriminator/DATA_LOCATIONS.md`
rows 18, 46-47) and DRAGEN (Timema, human, zebrafish benchmark runs; rows 11-17). For the manuscript, Timema
(curated v4 bank, `/staging/tmp/timema_sines/v4_curated/`) and *Scalopus* are the two that should be rerun with
the current code before any Results numbers are quoted; both need DRAGEN or KIT, or their genomes copied to
therioserver (Timema 1.24 GB on KIT, `/data/W/toki/Genomes/lower/Arthropoda/Timema/timema.fna`).
