# mme_main_20261024

*Mormoops megalophylla* (Mormoopidae; Noctilionoidea, Yangochiroptera), genome `GCA_985614085.1`. Complete SINEderella rerun on
main 1cd0b9b: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/mme_main_20261024/` (run 2026-10-24 14:28:26 to 2026-10-24 15:54:44; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 395 | 366 | 29 | 0 | 477 |
| VES | VES | 152,961 | 128,806 | 10,347 | 13,808 | 152797 |
| MEG-RL | MEG | 14 | 12 | 2 | 0 | 13 |
| MEG-RS | MEG | 160 | 70 | 5 | 85 | 173 |
| MEG-T2 | MEG | 49 | 49 | 0 | 0 | 56 |
| MEG-TR | MEG | 18 | 18 | 0 | 0 | 21 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 395 | 366 | 29 | 0 |
| VES | 152,961 | 128,806 | 10,347 | 13,808 |
| MEG | 241 | 149 | 7 | 85 |

Copies in the list 166,548; assigned 153,597; not SINE in total (assigned + not assigned) 26,844.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
