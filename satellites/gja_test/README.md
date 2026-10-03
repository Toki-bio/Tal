# Satellite screen, first real test: *Gekko japonicus* (KIT, 2026-10-03)

`tools/satellite_screen.py` (SINEderella `docs/SATELLITES.md`) on hits of the Squam3 consensuses, genome `Gekko_japonicus_V1.1`.

* `gja.*`: the existing hit BEDs of the user's earlier searches (Squam3A 216 438, Squam3B 185 230, sq3di 173 427 hits; 80 % length rule applied, so partial monomers are missing): kind A almost none (0.0-0.1 %), kind B 11-15 % of hits but 4-7 % expected by chance (excess 7.5-10 %): not flagged.
* `gja_partial.*`: Squam3A searched again with `sear Squam3A.q gja.bnk 0.2 65 0` (hits from 20 % of the consensus length, as in Vassetzky et al. 2023): 430 322 hits; **270 monomer runs (kind A, >= 4 monomers < 100 bp apart) holding 1 289 hits (0.3 %)**, median monomer (hit) length 141 bp (p10 112, p90 259), longest runs 56, 22, 10 monomers; 194 runs have 4 monomers. Kind B: 29.4 % in regular runs, 28.3 % expected by chance (excess 1.1 %): not flagged (the genome has one Squam3A hit per 6 kb).
* Reading: the SINE-derived satellite loci of the paper are found by coordinates alone, in seconds, but they are 0.3 % of the family's hits, so a family-share threshold can never flag Squam3A in this gecko; exclusion has to act on loci. Monomer type (which part of the SINE) needs the alignment coordinates, i.e. the characterisation step.
