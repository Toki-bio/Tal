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

1. **Layer 1: clear the major families first**: Dip, B2 (as written; see the note below), and perhaps B4. Each of the three is treated separately and he wants
   **a SubFam of 30,000 copies of each of the three families alone**; he then tells which subfamilies each consists of.
2. **Layer 2: deplete the genome of these majors**, treat the remaining genome the same way: re-run SINEderella to remove what is left of the three majors plus
   several next significant families.
3. Re-analysis with updated consensuses, then a re-evaluation of how the consensuses relate (overlaps).

*Note on names.* He wrote "dip, b2, and maybe b4". In this run the big rodent family is **B1** (538,486 firm), B2 has 1,717 firm / 10,294 soft and is
one he calls "OK old SINE, matches the original consensus poorly". To confirm with him whether the three majors are DIP, B1, B4 or DIP, B2, B4.

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

* **Extended flanks of top100 copies (B1, B1-dID, pB1).** The publish step's border loop extends a plate's flanks while the copies stay similar where the flank ends
  (`publish_sbe58.log`: "copies still similar where their flank ends - extending by +150/+0 bp" for several families). An extension can therefore reflect a shared flank
  (segmental duplication, tandem neighbours, a nested or composite element) rather than the SINE itself. The flank-twin table of this run (`results/flank_twins.tsv`;
  B1-dID 12.5 % twins, pB1 8.0 %, B4 9.3 %) is the first thing to check against the plates. **It has no row for B1 or DIP:** stage 9 (`fs9_twins.sh`) was killed
  by the system (gawk, out of memory) on the 538,486 B1 and the 1,185,121 DIP copies (`results/flank_twins/{B1,DIP}/fs9.log`; the run log says "B1 failed", "DIP failed"
  and the report simply has no row for them). The twin question for B1, the family he asks about, is therefore unanswered; it needs stage 9 in chunks (or on a
  sample of the family), which is a code change not made yet.
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
