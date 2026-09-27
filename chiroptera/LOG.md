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

## Hipposideros larvatus (hla) and Tadarida brasiliensis (tbr) — started 2026-09-26

His request: 1-2 bat genomes outside Vespertilionidae and Rhinolophidae, searched with SINEbase
Rhin-1 and VES; publish like the Tal species; if a SINE has >1000 copies, run de novo mode.

Genomes (chosen from his Hao et al. 2023 figure: Rhin hybridisation in Hipposideridae, Ves hybridisation+PCR in Molossidae):
| code | species | family | accession | level |
|---|---|---|---|---|
| hla | *Hipposideros larvatus* | Hipposideridae | GCA_031876335.1 mHipLar1.pri | chromosome, NCBI reference, contig N50 35 Mb |
| tbr | *Tadarida brasiliensis* | Molossidae | GCA_030848825.1 DD_mTadBra1_pri | chromosome, NCBI reference, contig N50 86 Mb |

Bank `~/chiro/chiro_bank.fa`: SINEbase Rhin-1 (182 bp) and VES (220 bp) with IUPAC codes resolved to a base
(his choice, 2026-09-26), so neither loses bases in sanitizing (earlier rsi/rre/rda/lly runs used a 177 bp Rhin-1).
- Rhin-1 Y,R,Y,R,W -> C,G,C,G,A: the bases of his rsi peel consensus_83, the group the Rhin-1 anchor fell in.
- VES R,Y -> G,T: SINEbase VES equals literature Ves_a (Fig. 14A) except those two codes; Ves_a has G and T there.

Run: therioserver `~/chiro/chiro_run.sh`, `SINEderella --publish`, THREADS=24 CHUNK_BP=30000 FLANK=50
SKIP_CANONICALIZE=1, both genomes in parallel.

"De novo mode" = his 2026-09-19 definition: AnnoSINE2 + SINE-de-novo-genome-scan, k-mer thinned, SubFam,
stopping at the candidate alignment for his review (`run_sinederella_combined.sh` de novo arm).

### Extended to one genome per family (his request 2026-09-26: "representative for each lineage")

Every family tip on his Hao et al. figure except Vespertilionidae and Rhinolophidae (excluded by the original request).
Choice rule: NCBI-flagged reference genome, highest contig N50, chromosome level preferred. hla and tbr above cover
Hipposideridae and Molossidae. Megadermatidae now uses *Macroderma gigas* (contig N50 38 Mb) instead of the
fragmented lly assembly. List with accessions: therioserver `~/chiro/chiro_genomes.tsv`.

| code | species | family | accession | assembly |
|---|---|---|---|---|
| rle | *Rousettus leschenaultii* | Pteropodidae | GCA_053538355.2 | ASM5353835v2 |
| tpe | *Triaenops persicus* | Rhinonycteridae | GCA_985607305.1 | mTriPer1.1.hap1 |
| mgi | *Macroderma gigas* | Megadermatidae | GCA_055471755.1 | mMacGig1.2_20250303 |
| cth | *Craseonycteris thonglongyai* | Craseonycteridae | GCA_985607335.1 | mCraTho2.1.hap1 |
| rmi | *Rhinopoma microphyllum* | Rhinopomatidae | GCA_043880545.1 | ASM4388054v1 |
| mau | *Myzopoda aurita* | Myzopodidae | GCA_985615965.1 | mMyzAur1.1.pri |
| nth | *Nycteris thebaica* | Nycteridae | GCA_985607365.1 | mNycThe1.1.hap1 |
| rna | *Rhynchonycteris naso* | Emballonuridae | GCA_037038555.1 | mRhyNas2.hap2 |
| tni | *Trinycteris nicefori* | Phyllostomidae | GCA_985607355.1 | mTriNic1.1.hap1 |
| mme | *Mormoops megalophylla* | Mormoopidae | GCA_985614085.1 | mMorMeg1.1.pri |
| nle | *Noctilio leporinus* | Noctilionidae | GCA_985617315.1 | mNocLep2.1.pri |
| fho | *Furipterus horrens* | Furipteridae | GCA_987537825.1 | mFurHor1.1.pri |
| ttr | *Thyroptera tricolor* | Thyropteridae | GCA_985616965.1 | mThyTri1.1.pri |
| mtu | *Mystacina tuberculata* | Mystacinidae | GCA_985612315.1 | mMysTub1.2.pri |
| cse | *Cistugo seabrae* | Cistugidae | GCA_061228675.1 | mCisSea1.1.hap1 |
| msc | *Miniopterus schreibersii* | Miniopteridae | GCA_964146895.2 | mMinSch1.hap1.2 |
| ntu | *Natalus tumidirostris* | Natalidae | GCA_985614685.1 | mNatTum1.1.pri |

Queue: `~/chiro/chiro_queue.sh` downloads all 17, waits for hla/tbr, then runs `run_one.sh` 3 at a time,
THREADS=20 each (60 total), same bank and settings as hla/tbr.

**Guard bug found:** `lly_run.sh` (and my `chiro_run.sh`, copied from it) called `disk_guard.sh 60 5`, but the first
argument is the PID to protect. PID 60 is a root kernel thread, so `kill -0` fails for toki and the guard exits
immediately: the lly run and the first minutes of hla/tbr ran with no disk guard. Started a correct guard on the
hla/tbr runner (PID 157826); the queue passes its own PID.

### Step1 hits, hla and tbr (2026-09-26, `genome.clean_step1/gen-consensuses.clean.part_<SINE>.bed`)

| code | Rhin-1 | VES | merged loci |
|---|---|---|---|
| hla | 70,888 | 55 | 70,910 |
| tbr | 452 | 641,631 | 640,752 |

Both exceed 1,000 copies -> de novo mode applies to both (`~/chiro/chiro_denovo.sh`, not yet launched: thread budget).

### De novo scope (his choice 2026-09-26, option 2)
De novo only for hla and tbr for now; he reviews those two candidate alignments before the other 17 get one.
Thread budget (cap 64): queue searches 3 x THREADS=10; each de novo chain AnnoSINE2 8 + scan 2 jobs x 4 = 16.
`~/chiro/denovo_launch.sh` starts both chains when the hla/tbr searches finish, each with its own disk guard.

### tbr published (2026-09-26)
run `~/chiro/tbr/run_20260926_223621`, publish rerun with RAW_ALN_BASE (MSA-viewer links).
| SINE | firm | soft | total | leak % | sim median |
|---|---|---|---|---|---|
| VES | 620,838 | 19,453 | 640,291 | 0.00 | 0.49 |
| Rhin-1 | 368 | 93 | 461 | 50.00 | 0.18 |
The discriminator calls both "SINE". Caution on Rhin-1: half its copies have VES within 10 % of the best score
(leak), and sim to Rhin-1 is 0.18, so these may be VES copies matched through the shared tRNA-derived head.
Not a verdict; to be judged on `tbr_Rhin-1_top100.aln.fa`.

hla: Rhin-1 70,868 (20,518 firm; 50,370 soft = rejected_low_bitscore, median copy 156 bp vs 182 bp consensus;
same 0.45 x 10th-best rule as rsi, 10th best 1,683 vs rsi 1,653). VES 42 copies, discriminator "No element".

### Page layout (2026-09-27)
Main page links to `chiroptera.html` (card + button beside Sicista); bat cards and status rows removed from
`index.html`. `chiroptera.html` is generated by `chiroptera/build_lineages.py` with the Tal page's own design:
its `<style>` block is copied verbatim from `index.html` at build time; sections are Cross-species resources,
one Species-style card grid per superfamily (tree order), and an Analysis Status table in the Tal format.

### Family tree on the page (2026-09-27)
Top card of `chiroptera.html`: Hao et al. 2023 Fig. 3 redrawn as inline SVG by `build_lineages.py` (TREE constant:
topology, node ages, crown ages of collapsed families; clade band and suborder bar colours matched to the figure; ICS
epoch colours). Each leaf has a Rhin (blue #0070C0) and VES (red #E00000) badge - the colours of the annotated figure:
solid >= 1,000 copies, pale 1-999, NA not searched, ... running; max over the family's genomes, hover lists all.
Caution recorded: badges encode copy number only. cth VES 96,273 copies at sim 0.18 and rmi VES 27,159 at sim 0.12
are solid red but far below tbr's 0.49; they need checking on the alignments before being read as VES presence.

### Vespertilionidae added (his request 2026-09-27)
Two distinct vespertilionid lineages, as positive controls for VES: *Myotis evotis* (Myotinae, mev,
GCA_056320405.1 mMyoEvo1.0, chromosome, NCBI reference, contig N50 101 Mb - the best vespertilionid assembly)
and *Vespertilio murinus* (Vespertilioninae, vmu, GCA_963924515.2 mVesMur1.hap1.2, chromosome, reference,
contig N50 49 Mb; the family's type genus). `~/chiro/vesp_queue.sh`: same bank and settings, starts when the
main queue ends (thread budget), 2 in parallel.

### 2026-09-27 fixes applied to all bat genomes
All 25 genomes are published. The 13 published before the fix were republished with the orientation and
shared-flank corrections (SINEderella ae764db, SINE-discriminator cedfe1b); later ones were published with them.
The only edge change in the republish: cth Rhin-1 rand100 +56 bp at 3'. Assembly-QC tables in the run dirs
were recomputed after a QC bug (SINEderella dd552e5); all 19 chromosome-level assemblies are SOLID.

### De novo candidates, hla and tbr (finished 2026-09-27 06:16 / 06:46 MSK)
`~/chiro/chiro_denovo.sh`: AnnoSINE2 mode 3 (`-a 2`, 8 threads, `-temd` on /home, HMM watchdog - no stalls)
+ sine_scan.sh (2a1dfa0, v2 SINEbase+literature fragment bank, 1764 fragments) -> merge -> k-mer thinning
(kmer_thin_singletons.py) -> SubFam 50. Stops at the candidate alignment for his review.
| code | AnnoSINE2 seeds | scan candidates | combined | kept after thinning | SubFam chunks |
|---|---|---|---|---|---|
| hla | 220 | 13,866 | 14,086 | 5,994 | 119 |
| tbr | 1,097 | 20,318 | 21,415 | 9,532 | 190 |
Published: `hla/alignments/hla_denovo_candidates_119chunks.aln.fa`, `tbr/alignments/tbr_denovo_candidates_190chunks.aln.fa`.
Clustered, not assigned or classified: input for manual review, not a result.

### MEG families added (2026-09-27)

Gogolevsky, Vassetzky & Kramerov 2009, Genomics 93:494–500. Four SINEs described from megabats and
reported absent from microbats and other mammals: MEG-RL (5S rRNA, paper consensus 213 bp), MEG-RS
(5S rRNA, 135 bp; nearly the 5′ part of MEG-RL, but copies carry their own TSDs), MEG-TR (tRNA/5S hybrid),
MEG-T2 (tRNA-derived, variable CNCCRGG tandem region). Paper copy estimates in *Pteropus vampyrus*:
about 10^4 MEG-RL, 7×10^3 MEG-RS, 5×10^3 MEG-TR, 4×10^3 MEG-T2.

Taken from the SINEbase download `~/SINEs_post2015.fas`, not retyped from the paper figures.
Lengths in that file: MEG-RL 207 bp (six shorter than the paper's 213; the SINEbase record is what was
searched), MEG-RS 135 (matches the paper), MEG-T2 232, MEG-TR 196. MEG-RL, MEG-RS and MEG-TR are ACGT.
MEG-T2 had three IUPAC codes, resolved the same way as Rhin-1 and VES (a base, not a deletion):

| site | code | base | reason |
|---|---|---|---|
| 121 | K (G/T) | T | the paper's inter-repeat linker is written `tg` (`CCCCRGGtgCGCCAGG`) |
| 127 | R (A/G) | A | that same linker continues `tgCGCCAGG`, so the R in `CGCCRGG` is A |
| 172 | R (A/G) | A | the paper's tail is AT-rich; A continues the A stretch |

Bank: `~/chiro/MEG.resolved.fa`. Names kept as MEG-RL, MEG-RS, MEG-T2, MEG-TR.
`SKIP_CANONICALIZE=1`: canonicalize merges seeds at 80% identity and rewrites every sequence's orientation.
MEG-RS is almost MEG-RL, and the old Rhin-1 and VES consensuses must stay byte-for-byte the ones already searched.

Add, not a new search. `SINEderella --add` copies the finished run, runs `sear` only for the new names
(0.8, 65, 50, same as step 1), rebuilds the merged intervals, and re-votes only loci that overlap a new
hit plus anything that was not a 10/10 assignment. Old 10/10 assignments that do not overlap a MEG hit
are carried forward. `~/chiro/meg_add.sh`, three genomes at a time, THREADS=10.
*Rousettus leschenaultii* (rle) is the positive control: the paper found MEG-RL in this species.
A microbat with a real VES signal (mev or tbr) is the negative control: the paper found no MEG outside Pteropodidae.

The page is still `build_lineages.py`. A MEG column and a purple tree badge appear only once a summary
contains the four names; until then the badge is NA, not a false zero. The badge is the sum of the four
assigned totals. The old SubFam plate is not rebuilt: `--add` does not rerun SubFam, so that button stays
the Rhin-1/VES plate.

### MEG results pushed (2026-09-27): rle, tbr, tpe, hla, mgi, mev, cth, rmi

Search of the four new consensuses is the fast part: on mev, `sear` for all four took 4 minutes and the
incremental re-vote 2 minutes. The wall clock after that is publish. On mev the border loop walked MEG-RS
5′ out to the 1000 bp cap and it was still shared (`STILL_BAD`, island fraction 0.33): 54 minutes. step8a
then aligned those widened copies (median 1640 bp around a 135 bp consensus) with
`mafft --localpair --maxiterate 1000`: 9 minutes for MEG-RS alone. Whole add: mev 79 min, tbr 43, hla 23,
rle 10, tpe 9, mgi 7. On hla the border loop raised `IndexError` in `border_loop_subfam.py` for two small MEG
sets; the pipeline logged WARN and published those plates without extension.
The queue is two adds at a time because `disk_guard` holds a job slot.

| code | MEG-RL | MEG-RS | MEG-T2 | MEG-TR | note |
|---|---|---|---|---|---|
| rle | 9,897 sim 0.73 | 7,383 sim 0.72 | 2,324 sim 0.41 | 2,775 sim 0.80 | positive control; order of magnitude matches Gogolevsky 2009. MEG-RL conflict 95%, MEG-TR 94%: MEG-RS is nearly the 5′ of MEG-RL, so the four counts are winning labels, not four clean populations |
| tbr | 14 sim 0.16 | 273 sim 0.79 | 92 sim 0.07 | 26 sim 0.16 | MEG-RS at 0.79 is not background. Look at `tbr_MEG-RS_top100.aln.fa` before calling MEG absent from this molossid |
| tpe | 16 sim 0.19 | 57 sim 0.24 | 118 sim 0.09 | 26 sim 0.12 | MEG here is weak |
| hla | 4,034 sim 0.27 | 113 sim 0.63 | 90 sim 0.07 | 766 sim 0.31 | MEG-RL count is large and the similarity is low; MEG-TR conflict 97% |
| mgi | 19 sim 0.18 | 310 sim 0.20 | 70 sim 0.09 | 39 sim 0.14 | background; MEG-RS firm only 9 of 310 |
| mev | 8 sim 0.21 | 372 sim 0.77 | 128 sim 0.11 | 27 sim 0.14 | MEG-RS like tbr: high similarity, and its 5′ flank is shared past 1000 bp, so these hits sit inside a longer common sequence. Not called a MEG family until the plate is read |
| cth | 21 sim 0.12 | 245 sim 0.19 | 68 sim 0.08 | 34 sim 0.14 | background; MEG-RS firm only 4 of 245 |
| rmi | 8 sim 0.13 | 162 sim 0.19 | 222 sim 0.09 | 87 sim 0.13 | background |

On the page, NA means that SINE has not been searched in that genome yet; a searched SINE with no copies shows 0.
MEG NA marks genomes still in the add queue. VES NA on Rhinolophidae is permanent for the current runs:
rsi, rre and rda were searched with Rhin-1 only. lly had no plates because its report was built without publish;
its MEG add runs with `--publish`, which writes them.
The tree's badge column was moved right (Rhin 992, VES 1056, MEG 1120) and the superfamily labels and suborder bars
moved past it, because the MEG badge sat on the superfamily text and the three *Rhinolophus* names reached the Rhin badge.

Reports and the new plates are in each species directory. Rhin-1 and VES counts are the carried-forward assignments (tbr VES 640,224; tpe Rhin-1 56,208; hla Rhin-1 70,867).


- 2026-09-27: rre (*Rhinolophus rex*) MEG add published (run_add_20260927_181142): Rhin-1 32,408 (firm 23,305, sim median 0.50); MEG-T2 114, MEG-TR 53, MEG-RS 37, MEG-RL 19, all at sim median 0.10-0.22. Card now links rre/report.html (dbe04cd). rda still pending.
- 2026-09-27: rda (*Rhinolophus darlingi*) MEG add published; card links rda/report.html. rsi MEG add published (run in ~/rhin/rsi/run_add_20260927_180847, from the peel10 bank; bundler fixed to look there): MEG-RS 1,586 (sim median 0.22), MEG-TR 170, MEG-T2 116, MEG-RL 22. Only the MEG plates were copied; hand-corrected r1-r10 plates kept.
