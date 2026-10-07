# Sicista betulina, 58-consensus run (sbe58): his observations on the plates, 2026-10-07

Source: his reading of the published plates (top 100, 100 random, SubFam) of `sicista/sbe58/`, run `run_20261006_143837`
(SINEderella 2fac340). Counts are firm / soft from `sicista/sbe58/summary.by_subfam.tsv`. The right-hand column is his wording, shortened
only where he repeated himself; his open questions are kept as questions. Nothing here is a result of mine; the leads under section 4 are
unverified.

**Overall impression (his):** 5-7 solid families of different age, heavily dominated by Dip and the rodent-type families. Needs re-analysis with
updated consensuses and a re-evaluation of the consensus interrelations (could be overlaps). Families not listed below are probably
misattributed scattered artifacts.

## 1. What he sees, per family

| consensus | firm / soft | his call | open points he named |
|---|---|---|---|
| **DIP** | 1,185,121 / 10,277 | good SINE family, 3 or 4 subfamilies | - |
| **B1** | 538,486 / 42,915 | good SINE family, maybe 1 or 2 subfamilies; needs checking by assignment | why top100 copies sometimes have an extended left flank; why random100 copies have TC motifs in the left flank (nested insertions? composite element?); needs a family-wide or genome-wide check |
| **B4** | 90,751 / 74,266 | good old SINE, maybe 2-3 not very clear subfamilies | - |
| **B1-dID** | 41,134 / 8,305 | looks like a legitimate SINE, 2-3 subfamilies possible but the SubFam data are unclear | same left-flank extension problem; right flank badly resolved in random100 |
| **pB1** | 20,643 / 28,675 | good SINE family, 1 or maybe 2 subfamilies (unstable insertion in the middle) | some right flanks look extended |
| **vic-1** | 16,482 / 524 | very old but legitimate SINE family | the bottom 3 sequences of SubFam need investigation |
| **STRIDM** | 13,464 / 1,772 | very old, can be a legitimate SINE | right flank needs clarification; many strange sequences at the bottom of SubFam |
| **Tu-II** | 5,829 / 484 | strange: top100 and random100 look like a very old SINE | SubFam sequences match only the middle of the consensus: why? |
| **DAS-I** | 5,760 / 13,887 | good, but a very old SINE | - |
| **IDL-Geo** | 3,673 / 2,444 | left flank OK | right flank very problematic; maybe a bipartite SINE with an unstable structure of parts |
| **RSINE1** | 3,288 / 7,466 | old; divergence is strange in the middle and at the right end; can be a SINE | bottom SubFam sequences need clarification |
| **B2** | 1,717 / 10,294 | OK, old SINE | matches the original B2 consensus poorly: consensus needs revision |
| **CAN** | 1,558 / 1,565 | big problematic unknown mess, probably pulled in by a long TCTC stretch | to be treated separately |
| **MEN** | 1,294 / 4,528 | maybe a SINE | why does SubFam look so terrible? |
| **Mon-1** | 1,227 / 26 | probably an OK ancient SINE | SubFam very discordant |
| **TUB** | 782 / 2,279 | OK, but may be mistaken for another SINE | SubFam plagued by a TC region |
| **Mar3** | 709 / 623 | maybe OK | differs from the original consensus: to verify |
| **ERI-1** | 553 / 5,033 | OK | the middle is a bit unstable |
| **MyrSINE** | 522 / 14,245 | alignment looks good | totally different from the original consensus: maybe another family |
| **TAL** | 312 / 22 | no element signal | - |
| **Mar1** | 307 / 17 | something very weak but detectable | to re-clarify |
| **CYN-III** | 246 / 4 | very weak signal, no defined left border, maybe not a SINE | - |

## 2. His plan for the peeling (layered)

1. **Layer 1: clear the major families first**: Dip, B1 and perhaps B4. Each of the three is treated separately and he wants
   **a SubFam of 30,000 copies of each of the three families alone**; he then tells which subfamilies each consists of.
2. **Layer 2: deplete the genome of these majors**, treat the remaining genome the same way: re-run SINEderella to remove what is left of the three majors plus
   several next significant families.
3. Re-analysis with updated consensuses, then a re-evaluation of how the consensuses relate (overlaps).

*Names (his answer 2026-10-07).* "b2" in his plan was a typo: **the three majors are DIP, B1 and B4.** B2 (1,717 firm / 10,294 soft) is one of the small old families.

## 3. Is the two-layer scheme (families first, subfamilies within families) implemented? (read 2026-10-07, SINEderella 2fac340)

**No.** `step2_asSINEment.sh` treats the bank as a flat list: every consensus (family, subfamily, length version or composite) is an equal competitor, ten
`ssearch36 -z 11` cycles, unanimity rule, per-subfamily threshold. There is no family label in the bank, in step 2 or in the outputs. What exists:

* the evaluation `SINEderella/docs/FAMILY_SUBFAMILY_ASSIGNMENT.md` (2026-10-05): family calls re-tallied from the kept votes on rsi, hs21 Alu and mm19 B1/B2 are
  99.4-100 % firm and library-size independent; of 1,994 vote failures only 7 are across families; not implemented. It lists four decisions only he can make
  (what the families are in the bank, where composites go, the layer-2 rule, whether family-firm / subfamily-unresolved copies count);
* **genome depletion is implemented**: `SINEderella --mask-bed BED` and `--mask-run RUN[:NAME,...]` write the loci of the named families (all subfamilies together)
  as N into a working copy of the genome with unchanged coordinates (`docs/MASKING.md`; `tests/test_mask.sh` exists, not run in the 2026-10-05 release check).
  Layer 2 of his plan can therefore be run today by hand with the existing flat assignment.

## 4. Leads for his open questions (not verified, no verdicts)

* **Left flanks that stick out of the alignment (B1 top100; B1-dID, pB1 not measured yet).** His reading: the flanks protrude from the alignment (the opposite of
  shared, extended flanks). Measured on `sbe58_B1_top100.aln.fa` (100 copies): 12 rows carry 79-96 bases to the left of the columns where at least half of the rows have
  sequence, the other 88 rows carry none (the plate columns 0-95 are occupied by 11-18 rows). These 12 copies are **longer loci**: 448-680 bp against a median of 364 bp
  for the other 88 (359-579). Their protruding left sequences are **unrelated to each other** (pairwise identity median 25 %, range 16-35 %, the background level), so they
  are not a shared flank, and they lie on 12 different contigs (7 of the 12 carry the `[array]` tag). The same measurement on `sbe58_B1_rand100.aln.fa` finds no row with
  more than 18 bases left of the core (6 rows locus > 470 bp on the right). So what sticks out is extra sequence belonging to the locus (nested insertion, composite element or a
  merged neighbour), not an artefact of the publish border loop. *My first lead (shared flanks from segmental duplication, extended by the border loop) was not supported
  and is withdrawn.* Still open: what that extra left sequence is (search it against the bank and the genome); the 12 rows are listed in the plate by name.
* **Flank-twin check (stage 9) for B1 and DIP.** The run's table has no B1 or DIP row: stage 9 was killed by the system (gawk out of memory) on the 538,486 B1 and the 1,185,121
  DIP copies. Rewritten (SINEderella d418439: the candidate step streams) and tested: identical output to the old script on the real B1-dID (40,941 copies), 11.9 GB -> 1.6 GB.
* **TC motifs in left flanks (B1 random100, TUB SubFam, CAN).** Possibly a low-complexity or microsatellite-prone neighbourhood; the satellite screen only removes kind-A
  satellites of the SINE and kind-B arrays, not simple repeats beside copies.
* **Consensus mismatch (B2, Mar3, MyrSINE).** The plates' row 1 is rebuilt from the copies and compared with the bank consensus in the consensus audit
  (`results/consensus_audit/`); its verdicts for these families are the starting point for a consensus revision.
* **SubFam sequences matching only the middle (Tu-II) or badly resolved (MEN, Mon-1, STRIDM, vic-1 bottoms):** SubFam works on the 30,000-copy sample drawn from
  all families together; per-family samples (point 1 of his plan) remove that mixing.

## 5. What a per-family 30,000-copy SubFam needs (proposed, not started)

The run kept every hit: `run_20261006_143837/results/assigned.fasta` (and `assignment_full.tsv`) hold the copies per family; the search is not repeated. For DIP, B1 (or B2) and B4:
select the family's firm and soft copies, draw 30,000 at random (seed 42, as for the other pages), SubFam (chunks of 50, 600 rows) and publish the alignment next to the
others. Cost: minutes to an hour per family on therioserver. Then his subfamily calls, then updated consensuses, then `--mask-run` of the three families and the
re-run for the remaining genome. Waiting for his go and the DIP/B1/B2/B4 naming above.
