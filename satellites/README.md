# satellites/ — SINE-derived satellites and SINE-containing tandem arrays

Tests and results of the SINEderella satellite screen (design and status: SINEderella `docs/SATELLITES.md`). Browser: [viewer.html](https://toki-bio.github.io/Tal/satellites/viewer.html)
(loads any `loci.bed` of a run's `results/satellites/` or a `*.loci.tsv` of `tools/satellite_trf_verify.py`; summary cards, where the monomers sit on the SINE, monomers per locus, sortable locus table).

| folder | what |
|---|---|
| `gja_test/` | first test on *Gekko japonicus* hits: the position-only screen (superseded for kind A) |
| `nna_test/` | Indian cobra, Squam3C: genome-wide TRF ground truth (2 580 sSat3 loci) vs the position-only screen |
| `gate_trf_test/` | cobra and gecko with the gate + windowed TRF stage (cobra 89.7 % of the ground truth in 29 s): `*_final.loci.tsv`, `*.units.fa` |
| `bats_test/` | stage on the tbr and rle runs (tables only): tbr MEG-RS = tandem array; rle MEG-RS dispersed with long arrays of the rsi 2 155 bp unit inside |

Paper: Vassetzky, Kosushkin & Ryskov 2023, Mob DNA 14:21, doi:10.1186/s13100-023-00309-2.
