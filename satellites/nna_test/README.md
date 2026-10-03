# Satellite screen, second real test: Indian cobra (*Naja naja*, Nana_v5 = GCA_009733165.1), Squam3C (KIT, 2026-10-03)

Ground truth from Tandem Repeats Finder (`trf 2 5 7 80 10 40 300`, whole genome in 44 parts, ~9 min on KIT): records with period 55-150 bp and
>= 4 copies, kept when the repeat unit aligns to Squam3C (ssearch36 over >= 45 bp, E <= 0.01): `sat3_trf_loci.tsv` (contig, start, end, period, copies, TRF % match).

* **2 580 SINE (Squam3C)-derived tandem loci**, 21 799 monomers, 1.59 Mb. Period 67 bp in 1 750 of them (66: 156, 68: 64; 133-134: 113 = dimers of the 67 bp unit): the sSat3 unit of snakes (positions 42-108 of Squam3C, 3' part of the head with box B) of Vassetzky et al. 2023. Copies per locus: 4-9: 2 096, 10-19: 388, 20-49: 79, 50-99: 13, 100-199: 3, >= 200: 1 (559 monomers on SOZL01001521.1). Mean TRF match 83 %.
* **Screen from `sear` hits (Squam3C, hits from 20 % of the consensus, 54 269 hits), `nna.*`:** 156 monomer runs (kind A, >= 4 hits < 100 bp apart), longest 20, 1.5 % of hits; kind B excess 2.0 %: nothing flagged. Only **186 of the 2 580 TRF loci (7 %)** overlap a screen A run, but 2 261 (88 %) overlap some `sear` hit. So the hits are there; the rule "4 consecutive hits" fails because `sear` merges neighbouring monomers into one hit (hit length 69-76 bp, spacing 134 bp = two monomers per hit) and most loci have only 4-9 monomers.
* Consequence for the design (docs/SATELLITES.md 5b): hit spacing alone is not a reliable kind-A detector; use hit clusters as a gate and verify with TRF on those windows only.
