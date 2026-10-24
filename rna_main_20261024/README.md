# rna_main_20261024

*Rhynchonycteris naso* (Emballonuridae; Emballonuroidea, Yangochiroptera), genome `GCA_037038555.1`. Complete SINEderella rerun on
main 0bca33e: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/rna_main_20261024/` (run 2026-10-24 09:05:04 to 2026-10-24 10:24:21; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 853 | 639 | 195 | 19 | 878 |
| VES | VES | 120,311 | 80,231 | 19,341 | 20,739 | 121082 |
| MEG-RL | MEG | 8 | 7 | 1 | 0 | 25 |
| MEG-RS | MEG | 149 | 120 | 24 | 5 | 11 |
| MEG-T2 | MEG | 113 | 107 | 6 | 0 | 114 |
| MEG-TR | MEG | 19 | 16 | 3 | 0 | 16 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 853 | 639 | 195 | 19 |
| VES | 120,311 | 80,231 | 19,341 | 20,739 |
| MEG | 289 | 250 | 34 | 5 |

Copies in the list 121,614; assigned 121,453; not SINE in total (assigned + not assigned) 20,924.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
