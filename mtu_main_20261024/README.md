# mtu_main_20261024

*Mystacina tuberculata* (Mystacinidae; Noctilionoidea, Yangochiroptera), genome `GCA_985612315.1`. Complete SINEderella rerun on
main 8d4961c: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/mtu_main_20261024/` (run 2026-10-24 17:29:07 to 2026-10-24 18:34:43; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 333 | 308 | 23 | 2 | 381 |
| VES | VES | 161,596 | 133,879 | 13,659 | 14,058 | 161273 |
| MEG-RL | MEG | 7 | 7 | 0 | 0 | 8 |
| MEG-RS | MEG | 60 | 28 | 4 | 28 | 72 |
| MEG-T2 | MEG | 100 | 92 | 2 | 6 | 107 |
| MEG-TR | MEG | 14 | 11 | 1 | 2 | 26 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 333 | 308 | 23 | 2 |
| VES | 161,596 | 133,879 | 13,659 | 14,058 |
| MEG | 181 | 138 | 7 | 36 |

Copies in the list 177,555; assigned 162,110; not SINE in total (assigned + not assigned) 29,541.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
