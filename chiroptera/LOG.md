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
