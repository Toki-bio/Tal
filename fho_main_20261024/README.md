# fho_main_20261024

*Furipterus horrens* (Furipteridae; Noctilionoidea, Yangochiroptera), genome `GCA_987537825.1`. Complete SINEderella rerun on
main 0bca33e: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/fho_main_20261024/` (run 2026-10-24 10:44:10 to 2026-10-24 12:04:41; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 286 | 253 | 21 | 12 | 341 |
| VES | VES | 140,399 | 111,209 | 12,312 | 16,878 | 140383 |
| MEG-RL | MEG | 10 | 9 | 1 | 0 | 7 |
| MEG-RS | MEG | 74 | 26 | 7 | 41 | 141 |
| MEG-T2 | MEG | 128 | 112 | 7 | 9 | 133 |
| MEG-TR | MEG | 16 | 16 | 0 | 0 | 17 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 286 | 253 | 21 | 12 |
| VES | 140,399 | 111,209 | 12,312 | 16,878 |
| MEG | 228 | 163 | 15 | 50 |

Copies in the list 153,090; assigned 140,913; not SINE in total (assigned + not assigned) 29,117.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
