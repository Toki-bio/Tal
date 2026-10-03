# Gate + windowed TRF (`tools/satellite_trf_verify.py`), real tests on KIT, 2026-10-03

Hits (from `sear`, 20 % of the consensus length) are only a gate: a window (hit +- 1 kb, merged) around every hit, TRF on the windows, a TRF record is a locus when period >= 55 bp,
copies >= 4 and its unit (doubled) aligns to the SINE consensus over >= 45 bp. Files: `*_final.loci.tsv` (one row per locus with the aligned part of the SINE), `*.units.fa` (monomer per locus), `*.summary.tsv`.

| Genome / consensus | hits | windows | time | SINE-derived tandem loci | monomers | period mode | SINE part (mode) |
|---|---|---|---|---|---|---|---|
| Indian cobra / Squam3C (reconstructed) | 54 269 | 45 363 windows, 103.5 Mb (6 % of the genome) | 29 s | **2 052** | 15 949 | 67 bp (1 565) | 30-110 / 40-110 (paper: 42-108) |
| *Gekko japonicus* / Squam3A | 430 322 | 42 252 windows, 112 Mb (window per hit would be 727 Mb > cap 600 Mb: fell back to clusters of >= 2 hits) | 21 s | **537** | 2 865 | 117 / 114 / 126 bp | 100-220 |

Cobra against the genome-wide TRF ground truth (2 580 loci, `satellites/nna_test/`): **2 315 recovered (89.7 %)**: 1 889 of 2 096 loci with 4-9 monomers, 344 of 388 with 10-19, 82 of 96 with >= 20. Earlier settings on the same data: 100 bp padding 59.5 % (2 hits per window) / 62.5 % (every hit); 1 kb padding with 2 hits per window 61.4 %. Padding and "a window around every hit" both matter.
Toy genome (planted 30- and 6-monomer satellites, a 20-unit 1.8 kb array, 200 dispersed copies, dimers, trimers, a duplicated region): both planted satellites found with the planted SINE part (131-250), nothing else; the long-unit array is found by `satellite_screen.py` (kind B).
Not checked: the paper's per-species counts; whether the gecko loci are the paper's three sSat3 variants.
