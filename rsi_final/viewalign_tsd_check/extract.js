// pulls the TSD finder functions out of script.js into a module, unchanged
const fs = require('fs');
const src = fs.readFileSync(process.argv[2] || 'C:/work/MSA-viewer-new/script.js', 'utf8').replace(/\r\n/g, '\n');
const names = ['_maxOf','_tsdBase','_collectUngappedBasesWithColumns','_residueNumbers','_tsdFlanks','_sliceCollectedWindow','_tsdMismatchCount','_tsdGapFraction','_isTsdReferenceRow','_largestNonGapBlockBounds','_pickTsdReferenceSeq','_seqsForTsdBoundaryScoring','_boundaryFromReferenceSeq','_columnConservationScores','_windowMean','_findSineBoundaryColumns','_tsdTrimPair','_findBestTsdInFlanks','_findBestTsdWide','_findTSD'];
let out = '';
for (const n of names) {
  const i = src.indexOf('\nfunction ' + n + '(') + 1;
  if (i < 1) throw new Error('missing ' + n);
  let k = src.indexOf('{', i), d = 0;
  for (; k < src.length; k++) { const c = src[k]; if (c === '{') d++; else if (c === '}') { d--; if (!d) break; } }
  out += src.slice(i, k + 1) + '\n';
}
out += 'module.exports = { ' + names.join(', ') + ' };\n';
fs.writeFileSync(process.argv[3] || 'tsdlib.js', out);
console.log('ok', out.length);
