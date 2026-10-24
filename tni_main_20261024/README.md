# tni_main_20261024

*Trinycteris nicefori* (Phyllostomidae; Noctilionoidea, Yangochiroptera), genome `GCA_985607355.1`. Complete SINEderella rerun on
main 0bca33e: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/tni_main_20261024/` (run 2026-10-24 10:46:48 to 2026-10-24 13:24:40; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 250 | 223 | 27 | 0 | 320 |
| VES | VES | 274,816 | 214,478 | 21,342 | 38,996 | 273912 |
| MEG-RL | MEG | 10 | 10 | 0 | 0 | 9 |
| MEG-RS | MEG | 78 | 65 | 13 | 0 | 273 |
| MEG-T2 | MEG | 104 | 102 | 2 | 0 | 113 |
| MEG-TR | MEG | 15 | 13 | 2 | 0 | 16 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 250 | 223 | 27 | 0 |
| VES | 274,816 | 214,478 | 21,342 | 38,996 |
| MEG | 207 | 190 | 17 | 0 |

Copies in the list 289,102; assigned 275,273; not SINE in total (assigned + not assigned) 52,825.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
