# Apple (*Malus domestica*) TRF motifs after the simple-repeat filter

Stage 02-04 of satpipe (Toki-bio/CLsat_workflow, tag-less main at clone time 2026-11-26) on GCF_042453785.1 (GDT2T_hap1),
TRF `2 5 7 80 10 50 2000`, motif period 80-300, at least 10 copies, monomerize on, simple-repeat filter (internal coverage > 50 bp dropped).
Run on DRAGEN (Biotech2024): `/staging/tmp/apple_minisat/run_full/motifs/mdo.fa.clean.fas` (259 of 301 TRF motifs).

- `mdo_trf_filtered.unaligned.fa`: the 259 motifs as written by the pipeline (arbitrary start; not yet rotated to a canonical phase).
- `mdo_trf_filtered.aln.fa`: the same sequences sorted by length, `mafft --localpair --maxiterate 1000`.

Placeholder: the run used a dummy family consensus (Darevskia CLsat1) only so the config validates; it can change the header
coordinates of the 7 long monomerized motifs, not any motif sequence.
