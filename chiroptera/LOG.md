# Chiroptera SINEs

Bat genomes, one card per species on `chiroptera.html`. Rhinolophidae
(*R. sinicus*, *R. rex*, *R. darlingi*) and Megadermatidae (*L. lyra*) are both
on that page. The Rhin-1 peel table stays in `rhin.html`.

## Lyroderma lyra (lly) — 2026-09-26

Assembly [GCA_004026885.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_004026885.1/) MegLyr_v1_BIUU.
Server: therioserver `~/lly/run_20260926_215400`. Bank: SINEbase Rhin-1 and VES only
(`SKIP_CANONICALIZE=1`, `THREADS=24`, `CHUNK_BP=30000`, `FLANK=50`). The Rhin peel
consensuses were not used.

Assembly QC did not finish: the stock script bubble-sorts every scaffold in awk, and this
assembly has 1,902,801 scaffolds (N50 96,489; 94.75% of scaffolds are shorter than 1 kb).
QC was stopped and Step 1 started. The verdict is diagnostic only.

| Consensus | Firm | Soft | sim_ratio median |
|---|---|---|---|
| Rhin-1 | 290 | 0 | 0.16 |
| VES | 11 | 0 | 0.14 |

301 merged loci. No leak, no conflict. Report: `lly/report.html`.
SubFam of the 301 copies: 6 chunk consensi plus the two anchors, `lly/lly_subfam_input.aln.fa`.
