# rle_main_20261024

*Rousettus leschenaultii* (Pteropodidae; Pteropodoidea, Yinpterochiroptera), genome `GCA_053538355.2`. Complete SINEderella rerun on
main 1cd0b9b: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/rle_main_20261024/` (run 2026-10-24 13:57:11 to 2026-10-24 14:28:25; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 7 | 7 | 0 | 0 | 9 |
| VES | VES | 14 | 9 | 2 | 3 | 14 |
| MEG-RL | MEG | 9,700 | 8,416 | 1,058 | 226 | 9710 |
| MEG-RS | MEG | 5,539 | 4,735 | 480 | 324 | 6718 |
| MEG-T2 | MEG | 2,154 | 1,921 | 173 | 60 | 2159 |
| MEG-TR | MEG | 2,710 | 2,457 | 173 | 80 | 2703 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 7 | 7 | 0 | 0 |
| VES | 14 | 9 | 2 | 3 |
| MEG | 20,103 | 17,529 | 1,884 | 690 |

Copies in the list 21,043; assigned 20,124; not SINE in total (assigned + not assigned) 1,612.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
