# Satellite stage on the bat runs (tables only, no exclusion), 2026-10-03

`tools/satellite_stage.py --no-exclude` on the existing SINEderella runs (`~/chiro/tbr`, `~/chiro/rle`, `run_add_20260927_143901`): kind A on the
length-filtered hits (no raw hits in those old runs), kind B on the full-length hits. 75 s for tbr (641 631 VES hits), 51 s for rle.

| Genome | Consensus | full hits | kind A loci | kind B runs | excess over chance | flag | note |
|---|---|---|---|---|---|---|---|
| tbr (*Taphozous*) | MEG-RS | 286 | 0 | 2 | 66.4 % | **SAT_B** | two arrays of 81 and 109 units, unit 926 / 912 bp (the ~910 bp period seen on the plates 2026-09-28) |
| tbr | VES | 641 631 | 13 (92 monomers) | 40 858 | 8.2 % (56.9 % observed, 48.7 % chance) | SAT_A | a hit per 3 kb: regular spacing by chance is the rule; 221 runs with >= 50 units against 0 expected |
| rle (*Rousettus*) | MEG-RS | 19 470 | 0 | 35 | 6.9 % | - | **the family stays dispersed**, but 26 of the 35 runs have >= 20 units (9 with >= 50), chance expects none at this density; one run of 36 units has the **2 145 bp unit of the rsi MEG-RS array** |
| rle | MEG-TR | 3 947 | 0 | 28 | 15.3 % | - | 14 runs of 20-49 units, chance expects none |
| rle | MEG-RL, MEG-T2 | 10 616 / 2 329 | 0 | 0 / 1 | 0 | - | |

Run-length calibration (observed kind-B runs by units vs the chance null of the same hits at genome-average density):
rle MEG-RS 5-9: 6, 10-19: 3, 20-49: 17, >= 50: 9, chance: none; tbr VES 5-9: 30 498 (chance 37 331), 10-19: 8 200 (6 320), 20-49: 1 939 (180), >= 50: 221 (0).

Consequence (docs/SATELLITES.md 5f): the family-share flag is right for tbr MEG-RS and wrong as the only rule for rle, where real arrays sit
inside a dispersed family. The stage now also excludes **long runs** on their own: the smallest run length that chance cannot explain
(rle MEG-RS: 5 units, tbr VES: 50), per consensus (`--exclude-b long`, the default).
