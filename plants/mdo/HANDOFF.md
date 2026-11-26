# Handoff: apple (mdo) SINE search, state of 2026-11-26

For a Claude session on another computer. This repository is public: no connection details, keys or Repbase
sequences belong in it. Server access is described in the owner's private Claude config (skill `server-connections`,
repository `claude-config`); pull that first, and ask the owner if it is not available.

## Read first, in this order

1. This file.
2. `peel/README.md`: the owner's peel, the comparison, and every run result.
3. The live pages: https://toki-bio.github.io/Tal/plants.html and
   https://toki-bio.github.io/Tal/plants/mdo/sd/report.html
4. Before any SINEderella work: `docs/AIM.md` in Toki-bio/SINEderella, and the skills `sinederella`,
   `think-before-work`.

## What the owner asked for, and where it stands

| request (owner, 2026-11-26) | state |
|---|---|
| search the apple genome for SINEs on DRAGEN with AnnoSINE2, SINEderella and Repbase consensuses from KIT | done |
| Plants subsection in this repository, with a ViewAlign consensus plate | done |
| record and match his peeling; run SINEderella with his consensuses plus Repbase SINE-1_Mad and SINE-3_Mad | done |
| remove SINE-1_Mad, keep his consensuses, build plates, publish | done |
| fix the report faults in SINEderella (Opus agent), re-run, republish | done |
| this handoff | this file |

**Current result** (run `run_excl_20261126_031836`): 18,887 copies; SINE 12,331, grey zone 1,715, not SINE 4,841.
Queries: the 7 peeled consensuses and Repbase SINE-3_Mad.

## Where things are

Server `Biotech2024` (DRAGEN), all under `/staging/tmp/apple_sines/`:

| what | path |
|---|---|
| genome: NCBI GCF_042453785.1 (GDT2T_hap1), 19 sequences, 651,332,481 bp, md5 checked | `genome/apple.fna` |
| Repbase query set, 84 land-plant SINEs (licensed, keep on servers; source on KIT `/data/V/toki/RepBase.bnk`) | `repbase_plant84.fa` |
| SINEderella clone, at bd2a54e; SINEplot at 23234c4; tool environment from `environment.yml` | `SINEderella/`, `SINEplot/`, `env/` |
| AnnoSINE2 2.0.9 environment (outside this folder, used unpatched) | `/staging/tmp/kmer_experiment/envs/annosine` |
| toy test, first 10 Mb of NC_091661.1 | `toy/` |
| AnnoSINE2 whole genome: 496 seeds in `out/Seed_SINE.fa`; launch script `run.sh` | `anno_full/` |
| seed plate files | `publish/` |
| SINEderella run, 9 queries (b5c8a4a) | `sd_peel9/run_20261126_011507` |
| first removal of SINE-1_Mad | `sd_peel9/run_excl_20261126_020256` |
| **current run** (bd2a54e), source of the published pages; script `sd_peel9/excl2.sh` | `sd_peel9/run_excl_20261126_031836` |

An older SINEderella checkout at `/staging/tmp/SINEderella` has local changes and is far behind GitHub: not used,
not touched.

## In this folder

| file | what it does |
|---|---|
| `make_seed_plate.sh` | server script that made the 496-seed plate, the seed table and `counts.txt` |
| `build_page.py` | builds `../../plants.html` from the files here and in `sd/` (the run id is written in it) |
| `strip_repbase.py` | copies a report and its alignments into `sd/` without Repbase rows or sequence |
| `peel/` | owner's consensuses, comparison scripts and tables, peel log |

To republish after a new run: copy `results/report.html`, `results/alignments/`, `results/summary.by_subfam.tsv`,
`results/continuation/continuation.tsv` and `step2/step2_output/loci_verdict.tsv` to a temporary folder; run
`strip_repbase.py <report> <alignments> <repbase.fa> SINE-3_Mad sd SINE-1_Mad`; remake `sd/verdict_counts.tsv`
(query, verdict, count from `loci_verdict.tsv` columns 2 and 4) and `sd/continuation_check.tsv`; update the run id in
`build_page.py`; run it; check every file for Repbase entries before committing. `index.html` card numbers are edited
by hand. `publish_run.sh` needs `RAW_ALN_BASE` ending in `/plants/mdo/sd/alignments/`, `SPECIES_CODE=mdo`,
`SINEPLOT_PY`, and `env/bin` first on `PATH`.

Copies of SINE-5_Mad-like consensuses in the published alignments are expected: they are built from apple copies
(SINE-5_Mad was never a query).

## Open points

- `consensus_3seqs` is doubtful: 339 of 400 tested copies continue 5' into the same consensus, 303 copies are grey
  zone; the report table still calls it SINE. The owner has not decided.
- The owner was asked for the member lists of his 7 groups, to score the automatic peel methods fairly. No answer yet;
  group membership in `peel/` is rebuilt by similarity.
- One haplotype only. The owner's general rule is to use every assembly of the individual; a second haplotype of this
  assembly was not found in NCBI's chromosome-level list.
- Run checks reported SHORTER or DIVERGED for some consensuses (see the report's consensus audit); not acted on.
- SINEderella, not fixed: a satellite locus shared by an excluded and a remaining consensus is dropped, not
  reassigned, after `--exclude`; composite search pieces may carry over (not traced); `tools/composite_scan.py` and
  `tools/coseg_compare/prep_labeled.py` still shorten names the old way.

## Working rules that applied in this session

Close-cooperation mode: small steps, show the result, ask at each choice of method, data or scope. A local hook
allowed one pipeline launch per owner message and demanded a written session brief first. Launch server jobs detached,
record the process ID, stop only by that process tree.
