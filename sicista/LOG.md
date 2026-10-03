# Sicista — SINEderella database arm

## 2026-09-26 — step 1, Mammalia-only bank

- Genome: therioserver `~/Sicista2026/assembly/consensus.fasta` (Medaka assembly).
- Bank: the 58 SINEbase consensuses classed Mammalia (`alignments/sicista_consensuses.fa`, therioserver
  `~/sine_runs/Sicista/db_based/mammal_consensus.fa`). The first attempt with the full combined bank was too slow and was
  stopped (partial kept as `partial_fullbank_run_20260926_143745`).
- Run: `~/sine_runs/Sicista/db_based/run_20260926_150828`, THREADS=32, CHUNK_BP=30000, FLANK=50; stopped on purpose
  when step 2 started (15:08 -> 16:25).
- Result: 2,170,177 merged hits. `alignments/sicista_subfam_input_30k.aln.fa` = SubFam `input.clw.al`:
  600 chunk consensi of a 30,000-copy sample + the 58 bank consensuses (658 rows).
- Per-consensus raw hit counts are on sicista.html (overlapping; not copy numbers).

## 2026-10-03 — de novo scan, S. betulina (primary) and S. trizona

- Genomes (therioserver): *S. betulina* hifiasm **primary** `~/Sicista2026/assembly/primary_eval/hifiasm_primary.p_ctg.fa`
  (1,144 contigs, 2,992 Mb; N50 18.25 Mb, BUSCO C 99.2%), not hap1 and not the Medaka assembly used on 2026-09-26;
  *S. trizona* GCA_982266845.1 (mSicTri.1, 162 sequences, 2,798 Mb). Run directories `~/sine_runs/Sbet`, `~/sine_runs/Stri`.
- **Subsample:** scan and AnnoSINE2 were run on a seeded random 10% of each genome (300 x 1 Mb windows, `genome.sub.fa`, 300 Mb).
  The full-genome `sine_scan.sh` was started and stopped after 30 min: with the rodent genome the fragment search is hit-dense
  (a 5 Mb chunk took 17 min at 4 threads), projected >20 h per genome. **The full genome was not scanned.**
- Scan: `SINE-de-novo-genome-scan/sine_scan.sh`, query bank `step0bank.query.fa` v2 (1,782 fragments, SINEbase + literature),
  MIN_ID 65, MIN_COV 0.90, FLANK 50, 10 Mb chunks. Candidates: betulina 279,697, trizona 273,316.
- AnnoSINE2 (mode 3, `-a 2`, `-numa 5000`): seeds betulina 1,090, trizona 818 (about 2 h each).
- Merge scan + seeds, `kmer_thin_singletons.py` (betulina 280,787 -> 269,166; trizona 274,134 -> 262,623), then
  30,000 random copies (seed 42) plus all kept AnnoSINE seeds -> SubFam (BnkSz 50). The 30k sample is the standard one on these pages;
  without it the thinned sets would give ~5,000 chunks. Alignments: `alignments/sbet_denovo_30k.aln.fa` (621 rows),
  `alignments/str_denovo_30k.aln.fa` (616 rows), row names `<sp>_dn_NNN`; seeds in `*_annosine_seeds.fa`.
- Status: preliminary, pending manual review. Not a finished result.

## 2026-10-03 — de novo scan on the Dip-depleted genome (minus bank), whole genome

- Dip loci = every hit of the Dip bank (SINEbase DIP + literature Dip_a1, Dip_a2, Dip_b1, Dip_b2) from the full SINEderella
  run on each genome (`~/sine_runs/{Sbet,Stri}/dip/run_*/results/all_hits.labeled.bed`), merged strand-agnostically:
  betulina 1,221,306 loci / 258.1 Mb, trizona 1,244,136 loci / 263.5 Mb. Complement of the genome, pieces >= 100 bp:
  betulina 1,052,242 pieces / 2,729 Mb, trizona 1,076,261 pieces / 2,530 Mb (`minus.fa`, headers `chrom:start-end()`).
- `sine_scan.sh` on the whole minus bank (not subsampled; 25 Mb chunks, TARGET_GENOME = original genome for coordinates),
  same query bank and thresholds as the subsample scan (MIN_ID 65, MIN_COV 0.90). It took about 2 h per genome, against a projected
  >20 h for the full genome with Dip included. Candidates: betulina 92,161, trizona 98,540 (against 279,697 / 273,316 on the 10% subsample with Dip).
- Thinning (singletons dropped): 71,856 / 78,108 kept. 30,000 random copies (seed 42), SubFam (BnkSz 50). AnnoSINE2 not run on the minus bank.
  Alignments: `alignments/sbe_minus_denovo_30k.aln.fa`, `alignments/str_minus_denovo_30k.aln.fa` (600 rows each, `<sp>_mn_NNN`).
- Status: preliminary, pending manual review. Not a finished result.

## 2026-10-03 — Dip SINE full SINEderella runs published as species reports (sbe, str); page rebuilt

- Whole-genome SINEderella steps 1-4 + report on each genome with the five Dip consensuses (SINEbase DIP + literature Dip_a1, Dip_a2, Dip_b1, Dip_b2;
  as searched, after SINEderella removed IUPAC codes: DIP 206 bp, Dip_a1 224, Dip_a2 219, Dip_b1 240, Dip_b2 266; identical bank in both runs, `consensuses.clean.fa`).
  therioserver `~/sine_runs/Sbet/dip/run_20261003_044818` (S. betulina primary assembly) and `~/sine_runs/Stri/dip/run_20261003_044700` (S. trizona GCA_982266845.1).
- Loci entering assignment: 1,221,306 (sbe), 1,244,136 (str). Firm assigned: sbe 1,032,760 (Dip_a1 846,431, DIP 112,492, Dip_a2 73,613, Dip_b2 131, Dip_b1 93),
  str 1,036,280 (Dip_a1 845,550, DIP 117,228, Dip_a2 73,114, Dip_b1 211, Dip_b2 177). Dip_b1 and Dip_b2 are almost all soft assignments.
- Published with SINEderella `publish_run.sh` (step 7 boundary refinement, step 8a plates with 50 bp left / 70 bp right flanks, SINE-discriminator verdict columns, seed as row 2);
  files in `sbe/` and `str/`: `report.html`, `alignments/` (top 100, 100 random, SubFam for subfamilies with >= 400 copies), `summary.by_subfam.tsv`, `assignment_stats.tsv`, `consensuses.clean.fa`.
- `sicista.html` is now generated by `sicista/build_page.py` from these files (`sicista/_old_sections.html` holds the 2026-09-26 sections kept verbatim).
