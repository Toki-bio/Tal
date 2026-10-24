# msc_main_20261024

*Miniopterus schreibersii* (Miniopteridae; Vespertilionoidea, Yangochiroptera), genome `GCA_964146895.2`. Complete SINEderella rerun on
main 8d4961c: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/msc_main_20261024/` (run 2026-10-24 16:07:35 to 2026-10-24 17:50:13; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 377 | 347 | 30 | 0 | 428 |
| VES | VES | 248,134 | 193,539 | 26,478 | 28,117 | 248001 |
| MEG-RL | MEG | 11 | 11 | 0 | 0 | 12 |
| MEG-RS | MEG | 28 | 19 | 1 | 8 | 51 |
| MEG-T2 | MEG | 93 | 90 | 3 | 0 | 103 |
| MEG-TR | MEG | 16 | 16 | 0 | 0 | 17 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 377 | 347 | 30 | 0 |
| VES | 248,134 | 193,539 | 26,478 | 28,117 |
| MEG | 148 | 136 | 4 | 8 |

Copies in the list 263,504; assigned 248,659; not SINE in total (assigned + not assigned) 42,970.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
