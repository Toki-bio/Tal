# cth_main_20261024

*Craseonycteris thonglongyai* (Craseonycteridae; Rhinolophoidea, Yinpterochiroptera), genome `GCA_985607335.1`. Complete SINEderella rerun on
main 0bca33e: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/cth_main_20261024/` (run 2026-10-24 09:30:54 to 2026-10-24 10:46:47; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 181 | 25 | 154 | 2 | 197 |
| VES | VES | 91,807 | 52,316 | 23,986 | 15,505 | 96272 |
| MEG-RL | MEG | 23 | 19 | 4 | 0 | 8 |
| MEG-RS | MEG | 41 | 28 | 8 | 5 | 4 |
| MEG-T2 | MEG | 62 | 54 | 4 | 4 | 67 |
| MEG-TR | MEG | 24 | 17 | 7 | 0 | 27 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 181 | 25 | 154 | 2 |
| VES | 91,807 | 52,316 | 23,986 | 15,505 |
| MEG | 150 | 118 | 23 | 9 |

Copies in the list 92,186; assigned 92,138; not SINE in total (assigned + not assigned) 15,564.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
