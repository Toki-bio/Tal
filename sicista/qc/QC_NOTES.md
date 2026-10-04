# S. betulina (ZBS2026-Sb1) assembly QC and Illumina decision  (started 2026-10-03)

Server: therioserver ~/Sicista2026/qc/ (BAM, calls, models, sine). Assembly = hifiasm PRIMARY `assembly/primary_eval/hifiasm_primary.p_ctg.fa`.

## A. Assembly QC (all measured 2026-10-03 unless noted)
| metric | value |
|---|---|
| contigs / total | 1,144 / 2.992 Gb (no gaps; hifiasm contigs) |
| N50 / L50 | 18.25 Mb / 48 |
| N75 / N90 | 8.03 Mb / 1.93 Mb |
| largest contig | 84.3 Mb |
| contigs >=10 Mb / >=1 Mb / >=100 kb | 91 (2.11 Gb, 70.5%) / 292 (2.80 Gb, 93.7%) / 750 (2.98 Gb, 99.5%) |
| contigs <100 kb | 394 (14.3 Mb) |
| GC | 42.5% (contigs >=1 Mb: 36.5-54.3, median 42.1) |
| BUSCO mammalia_odb10 (n=9226) | C 99.2% [S 97.1, D 2.1], F 0.2, M 0.6 (S. trizona ref: C 98.3%) - audit 2026-09-29 |
| Merqury (ONT k-mers, relative only) | QV 58.3, completeness 93.4% |
| telomeres (TTAGGG, >=60 bp in terminal 10 kb) | 0 of 620 contigs >=200 kb have both ends; 30 have one end |
| GenomeScope2 (ONT 21-mers, error 1.2%) | haploid 2.69 Gb, het 0.50-0.57% (unreliable: ONT k-mers) vs assembly 2.99 Gb (+11%) |
| mitochondrial contamination | 83 mt-like contigs (1.37 Mb) must be removed from the nuclear set; 305 nuclear NUMT loci (see SICISTA_MITO_FINDINGS.md) |
| read mapping (minimap2 map-ont, all reads -> primary, 97 GB BAM, `qc/sbet_ont.bam`) | median depth 35.1x (10-kb windows; 5-95%: 15-46x); 75.5% of windows 0.75-1.25x of median |
| half-depth sequence | 10.8% of bp at MAPQ>=10, but 6% (~179 Mb in 494 contigs) is still half-depth when MAPQ0 reads are counted; the rest is low-MAPQ repeat |
| half-depth contigs = hemizygous X/Y | 81 contigs >=1 Mb at half depth have NO >=90%-identity counterpart in the assembly (so not unpurged haplotigs); 89.3% of their alignment to S. trizona goes to ONE chromosome (OZ418353.1, 148 Mb; 165x enriched) and 3.4% to OZ418354.1 (20 Mb) => sample is almost certainly male, X (+Y) ~175 Mb; exclude from autosomal PSMC |
| diploid-normal contigs | 292 contigs, 2,671 Mb (89%); 212 contigs >=1 Mb = 2.64 Gb. Largest: ptg000057l 84.3, 018l 75.0, 003l 62.5, 566l 58.8, 116l 56.8 Mb |
| reads | 13.5 M reads, 104.25 Gb (~35x of 2.99 Gb), N50 12.05 kb, mean Q15.8, >Q20 only 21.5% of reads; read headers non-standard (no basecaller model tags) |

## B. ONT-only variant calling tests (Clair3 v1.x, model r1041_e82_400bps_sup_v500; test set = 5 autosomal contigs, 337 Mb, 35x; `qc/calls/`)
- Chemistry: header/identity (98.3% median read identity on mtDNA) consistent with R10.4.1; model choice not verified by the sequencing provider (no pod5/model tag).
- Full 35x: 627,351 PASS het SNPs = 1,860 / Mb (0.19%); 5,160 hom-alt SNPs. 83% of het calls have AF 0.35-0.65 (median 0.46). Ti/Tv of PASS het SNPs = 1.64 (expected ~2.0-2.2 for real mammal SNPs -> some transversion-rich false calls, order of 15%).
- Disjoint read halves (each 17.5x, `sub_A`/`sub_B`): 659k het SNPs each; 87.8% of A's calls are in B (and vice versa) without filtering; 99.3% of the A∩B set is recovered by the 35x call set; 99.6% of 35x calls are in A∪B.
  Calls private to one half (80k each) have Ti/Tv 2.9 (transition-rich, typical of ONT 5mC-related C>T errors, or real hets sampled once) -> halves overestimate theta by ~17% vs the 35x set. Hard QUAL/AF filters change het density 7-fold (1,950 -> 261 / Mb) and are depth-dependent: do not use fixed QUAL thresholds across depths.
- PSMC (psmc -N25 -t15 -r5 -p 4+25*2+4+6; subset only = 11% of genome; mu 5.4e-9/gen, g=1 yr illustrative): the two independent halves (A, B) give nearly identical curves for t > ~30 ky (Ne ~90-200k, plateau ~150k through 0.1-2 My) and agree on a 10-30 ky trough (~20-30k); the 35x call set agrees for t > 50 ky but differs 2-5x for t < 40 ky (its theta is 15% lower). Plot: C:\Temp\sm\psmc\psmc_test2.png. => deep-time history robust; recent (<40 ky) history call-set dependent.
- SINE-like indels (Sniffles2 2.8.1 on the same subset): 7,060 PASS INS/DEL >=50 bp; 4,664 at 100-700 bp, of which 2,741 (58.8%) match the Dip SINE bank (Dip_a2 1,571; Dip_a1 894; DIP 273) = 8.1 per Mb (~24k genome-wide), 1,493 deletions + 1,248 insertions relative to the assembly, ALL heterozygous (expected: same individual). Reproducibility: ~89-90% of 100-700 bp calls from one 17x half are found in the other half; 98% of half-calls are in the 35x set.

## C. More tests (2026-10-04)
- **Hemizygous X as built-in false-positive control** (male: X het SNPs should not exist): Clair3 on the 53 half-depth contigs >=1 Mb (108.9 Mb, ~17x). 5 of them (14 Mb: ptg000139l, 160l, 198l, 253l, 317l) have >600 het/Mb = they behave as diploid/PAR/duplicated and are excluded. On the other 94.8 Mb: 163 het SNP/Mb unfiltered = 8.4% of the autosomal 17x rate (1,950/Mb); with AF 0.35-0.65: 14/Mb = ~1% of the autosomal AF-filtered rate. False hets had AF median 0.29 (true ~0.46). => ONT R10.4.1 + Clair3 sup model gives ~1% false hets once an allele-fraction filter is applied (upper bound; includes mismapped paralogs).
- **Assembly base accuracy, ONT-only estimate:** homozygous-alt calls (reads disagree with assembly) on the 337 Mb test set: 5,160 SNPs + 11,505 indels = 49/Mb = 4.9e-5 => QV ~43 (upper bound on error; Merqury QV 58 is relative only). Indels dominate (homopolymers) -> frameshift risk in gene models; self-polishing with the Clair3 1/1 calls is possible without Illumina.
- **Phasing (whatshap, ONT reads, test set):** 99.1% of het SNPs phased; phase blocks span most of each contig (e.g. ptg000566l 58.8 Mb in 2 blocks, NG50 48.8 Mb). Haplotype-resolved calls are feasible from ONT alone.
- **SINE-like heterozygous indels, whole assembly (Sniffles2, autosomal-class contigs 2.68 Gb):** 70,853 PASS INS/DEL >=50 bp; 46,423 at 100-700 bp with sequence; 25,313 (54.5%) match the Dip SINE bank (Dip_a2 14,437; Dip_a1 8,255; DIP 2,565) = 13,820 deletions + 11,493 insertions vs the assembly = ~9.4 per Mb, median length 244 bp. (Script label "per Mb 75" in the log is wrong - used test-set length.) These are within-individual heterozygous indels; only a population sample can show allele frequencies.

## D. Final genome-wide ONT-only PSMC (2026-10-04; `qc/final/`, plot `sicista_psmc_ont_only_2026-10-04.png`)
- Input: Clair3 (sup_v500) on 351 diploid-depth contigs (2.684 Gb; mt-like and X-class excluded), PASS het SNPs; callable mask = MAPQ>=10, depth 18-62 -> 2.589 Gb callable; 100-bp bins; het-bin fraction 13.5% (AF 0.3-0.7) / 14.5% (PASS only); theta/bin 0.103 / 0.120.
- psmc -N25 -t15 -r5 -p "4+25*2+4+6"; 20 bootstraps (AF-filtered set). mu 5.4e-9/gen and g=1 yr are ILLUSTRATIVE (not verified for Sicista).
- Ne (AF-filtered; bootstrap 2.5-97.5%): 2 My ~190k (177-211k); 1 My ~148k (144-157k); 300 ky ~166k (160-168k); 100 ky ~117k (109-126k); 50 ky ~71k (69-85k); 20 ky ~80k (45-96k); 10 ky ~12k (9-12k). Present-day N0 ~48k (AF) / 55k (PASS only).
- Shape: stable plateau 150-190k from >2 My to ~100 ky, decline to ~70k at 50-60 ky, a bump to ~130k at 25-40 ky, then collapse to ~10-30k by 15-20 ky ago. PASS-only vs AF-filtered curves agree for t > 60 ky; they differ in timing (up to ~10 ky) and amplitude of the last 40 ky (and PSMC has little power for the last ~10 ky). Absolute Ne scales with theta, which shifts ~15% between filters.
- Genome-wide result vs 337-Mb subset tests: same deep-time Ne (~150k), so subset tests were representative.

## E. Verdict on additional Illumina (see message for the short version)
- Not needed on this individual for PSMC, gene phylogeny, or within-individual SINE indels: ONT-only false-het rate ~1% after AF filter (hemizygous-X control), independent read halves reproduce PSMC and 88-90% of calls, bootstrap CIs are tight, BUSCO 99.2%.
- Needed for population-level goals regardless of platform: MSMC/MSMC2 beyond one genome (>=2-4 individuals), SINE insertion polymorphism frequencies, population genetics. Cheapest route: 15-20x PCR-free Illumina on 4-8 additional individuals against this assembly.
- Optional on this individual: 30x Illumina (~80-90 Gb) as orthogonal truth (QV, FP/FN of hets, homopolymer indel polishing). Self-estimated assembly QV ~43 from ONT hom-alt calls; indels dominate.
- Not verified: sequencing chemistry/model, mu and generation time for Sicista, sex (inferred male from X depth), recall of true hets.
