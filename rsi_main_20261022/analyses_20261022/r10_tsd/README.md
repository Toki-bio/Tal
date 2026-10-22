# r10_19seqs copies with their target-site duplications (rsi_main_20261022)

`rsi_r10_tsd.aln.fa` + `rsi_r10_tsd.bed` (ViewAlign track, alignment columns). Built by `tsd_aln.py` on therioserver
`/home/toki/tmp/tsd_r10/` from the 315 assigned r10 copies of `rsi_main_20261022` (padded 60 bp 5', 1 kb 3';
`tsdscan.py`).

- Copies: 5'- and 3'-complete r10 match (260), nearest TSD starting at most 250 bp past r10's end (48); copies with
  verdict not_SINE are left out (none here), grey ones are marked [grey]. Rows sorted by where the 3' TSD starts.
- TSD: the 5' flank ending at most 4 bp before the element, repeated on the 3' side; 10-20 bp, <= 10 % mismatches; the
  longest at the nearest position. Null: the 5' side of one copy against the 3' side of the next copy gives a TSD in the
  same window in 15 of 260 pairs, against 48 real. So about a third of the TSDs shown can be chance; TSD is a supporting
  fact, not proof.
- Columns: 1-40 left flank (lower case) | 41-60 5' TSD (UPPER, right-justified) | 61-64 bases between the TSD and the
  element (lower) | 65-662 element: r10 match and everything up to the 3' TSD (UPPER, MAFFT --localpair; r10's
  consensus spans 65-240) | 663-682 3' TSD (UPPER, left-justified) | right flank (lower). The gaps before the 3' TSD
  are there to put all TSDs in the same columns.

`rsi_r10_tsd.snapshot.json`: the same alignment as a ViewAlign snapshot with every 5' and 3' TSD base marked (yellow) and the
BED track. Open: https://toki-bio.github.io/MSA-viewer/?snapshotFile=https%3A%2F%2Fraw.githubusercontent.com%2FToki-bio%2FTal%2F80c8de9%2Frsi_main_20261022%2Fanalyses_20261022%2Fr10_tsd%2Frsi_r10_tsd.snapshot.json

## Is what follows r10 shared with other subfamilies? (`rsi_r10_vs_P18.aln.fa` + `.bed`, `post.py`, `post2.py`, `ctx.py`)
The 100 bp right after r10's end in 271 r10 copies, searched (ssearch36, >= 40 bp, E <= 1e-3) against up to 300 copies of
every other subfamily, each read 400 bp past its own end; control = 100 bp of the same copies' 5' flank:
- 2,644 copy pairs hit INSIDE r10_r8_P18 copies' elements, starting mostly at P18 positions 100-260 and ending at 180-340;
  hits past the ends of other subfamilies' copies are at the control's level (genomic repeats in their flanks);
- P18's middle piece (146-291) is a tandem array with a 10-bp period (self-identity 0.61 at 10, 0.68 at 20 bp): units like
  GTGAGCTGCC / GTGGGCTGCT; r10's consensus ends with the start of it (...GCTGCCGTGAGCAGCC).
So the r10 copies carry a variable number of these units past r10's end, then a G-rich / A-rich stretch, then the TSD.
The alignment: r10 and P18 consensuses, the 48 r10 copies with a TSD (element up to their 3' TSD) and 40 P18 copies (SINE
verdict, best bitscore), MAFFT --localpair; the track marks P18's three parts on the P18 consensus row.
