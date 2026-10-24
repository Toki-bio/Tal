# nle_main_20261024

*Noctilio leporinus* (Noctilionidae; Noctilionoidea, Yangochiroptera), genome `GCA_985617315.1`. Complete SINEderella rerun on
main 1cd0b9b: bank of the six SINEbase consensuses used for the published 2026-09-27 run (Rhin-1, VES, MEG-RL, MEG-RS,
MEG-T2, MEG-TR), with a family table (Rhin, VES, MEG). New against that run: a verdict for every copy from its flanks
(SINE locus / grey zone / not a SINE), plates drawn from SINE loci with packed flanks, the continuation check, families.

- Server: therioserver `/home/toki/chiro/nle_main_20261024/` (run 2026-10-24 12:58:13 to 2026-10-24 17:48:49; scripts `/home/toki/chiro/rerun_20261024/`).
- `loci_verdict.tsv.gz`: every copy, its subfamily, assignment status, verdict and reason (gzipped; plain text inside).

| subfamily | family | assigned | SINE | grey | not SINE | firmly assigned in the 2026-09-27 run |
|---|---|---|---|---|---|---|
| Rhin-1 | Rhin | 155 | 134 | 21 | 0 | 174 |
| VES | VES | 434,859 | 330,837 | 44,865 | 59,157 | 433613 |
| MEG-RL | MEG | 37 | 23 | 3 | 11 | 12 |
| MEG-RS | MEG | 270 | 79 | 28 | 163 | 249 |
| MEG-T2 | MEG | 97 | 95 | 0 | 2 | 101 |
| MEG-TR | MEG | 52 | 26 | 2 | 24 | 54 |

| family | assigned | SINE | grey | not SINE |
|---|---|---|---|---|
| Rhin | 155 | 134 | 21 | 0 |
| VES | 434,859 | 330,837 | 44,865 | 59,157 |
| MEG | 456 | 223 | 33 | 200 |

Copies in the list 459,086; assigned 435,470; not SINE in total (assigned + not assigned) 82,973.
Rules: SINEderella `docs/VERDICTS.md`.

Contents: `report.html` (full report), `alignments/` (top 100, random 100, subfamily plates), `flankscan/`, `satellites/`,
`length_variants/`, `consensus_audit/`, `continuation/` (where present), per-subfamily tables, `consensuses.clean.fa`, `manifest.txt`,
`pipeline_stages.txt` (stage times).
