# r5_r5_P48 in rfe: two divergence peaks = two sequence variants

The r5_r5_P48 divergence histogram of `rfe_final_20261115` has two peaks (~9-10 % and ~16-17 %).

**Halves scored separately** (ssearch36 -z -1, each copy vs consensus pieces 1-136 and 137-308): low group (256 copies, < 13 %) head 9.1 %, second unit 11.0 %; high group (223 copies, >= 14 %) head 14.9 %, second unit 16.7 %; median copy length 318 vs 321 bp; the r5_r5 second unit fits better than the r6 (r5_r6_P26 138-376) or r3 (r5_r3_P34 141-339) units in both groups. So the high peak is not copies with a different second unit.

**Separation test** (`r5r5_lo50_hi50.aln.fa`): only copies with verdict SINE (unique flanks; 400 of 495). `lo_*` = the 50 lowest-divergence copies of the low group (5.6-8.6 %); `hi_*` = the 50 copies nearest the high-peak centre 16.9 % (16.2-17.6 %; the high group's own "best" copies cannot be chosen without its own consensus, and the copies nearest r5_r5 sit at the dip). Row names = group_divergence_locus; `CONS_r5_r5_P48` = the bank consensus. MAFFT L-INS-i `--reorder`. List: `r5r5_sel.tsv`.

**Result:** the guide-tree order puts 40 lo rows together and 46 hi rows together. 18 columns are diagnostic (majority base >= 60 % in each group and different), spread over the head (consensus positions 40, 48, 108, 115-117, 120, 135) and the second unit (147, 170, 175, 181, 197, 211, 213, 244, 268, 275). Control: 200 random 50/50 splits of the same 100 rows give at most 2 such columns (mean 0.04). In every diagnostic column the low group carries the consensus base.

**Reading:** the high peak is a second r5_r5-type variant whose own consensus differs from the bank's r5_r5_P48 at >= 18 positions along the whole length; its copies are scored against the wrong consensus, so their divergence is inflated. It is a candidate subfamily for the bank (owner's decision).
