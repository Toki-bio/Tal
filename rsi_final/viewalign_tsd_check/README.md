# ViewAlign TSD finder: misses on rsi plates (2026-10-29)

A colleague reported that ViewAlign's TSD finder misses obvious TSDs. These scripts measure that, using the
finder's own code, copied unchanged out of `script.js`.

```
node extract.js /path/to/MSA-viewer/script.js tsdlib.js     # the finder as a Node module
LIBS=./tsdlib.js node real.js       # the 18 plates in ../stage8_plates against SINEderella's TSD rule
TSDLIB=./tsdlib.js node sim.js      # planted-TSD alignments, 10 scenarios
TSDLIB=./tsdlib.js node fp.js       # copies with no TSD: how many get a call
TSDLIB=./tsdlib.js node oldchecks.js   # two fixtures from the ViewAlign regression suite
```

`real.js` compares against an exhaustive search with SINEderella's has_tsd rule: 10-20 bp, at most 10 %
mismatches, the 5' copy ending at most 4 bp before the element and the 3' copy starting 0-30 bp after it
(element = the REF/proposed columns). Its null (copy i's 5' side against copy i+1's 3' side) finds 59. The
CONTROL row runs the finder with each copy's 3' side swapped for the next copy's, so every call of 10 bp or
more there is chance.

| build | oracle TSDs found (of 941) | missed | a different one | CONTROL calls >= 10 bp (of 3049) |
|---|---|---|---|---|
| v300 | 635 (67 %) | 200 | 106 | 108 |
| v301 | 838 (89 %) | 58 | 45 | 238 |

v301 changes: the 3' copy is looked for from 6 bases before the element's end to 30 past it, with only strict
pairs (>= 10 bp, <= 10 % mismatches) away from the first 0-3 bases and an exact match on any part inside the
element. The automatic boundary is trimmed to a run of 3 conserved columns, and columns are scored by
identity among the copies present. Example of a v300 miss: plate r8_83seqs_S, copy w43714,
AGAAAAACGCCAGGTA / AAAAAAACGCCAGGTA after a long (AAATA)n tail.
