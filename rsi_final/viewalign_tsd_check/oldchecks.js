// the two existing regression fixtures (Repeat Finder planted 6-bp, flank check), auto mode, no browser
const T = require(process.env.TSDLIB);
function mk(seed0, bodyLen, n, plantedNoTsdFrom) {
  let seed = seed0; const rnd = () => Math.floor((seed = (seed * 1103515245 + 12345) % 2147483648) / 65536);
  const rb = k => Array.from({ length: k }, () => 'ACGT'[rnd() % 4]).join('');
  const body = rb(bodyLen); const rows = [{ header: 'consensus', seq: '-'.repeat(31) + body + 'AAAAAAAAAA' + '-'.repeat(31) }];
  for (let i = 0; i < n; i++) { const tsd = rb(6); const b = body.split('').map(c => (rnd() % 20 === 0 ? 'ACGT'[rnd() % 4] : c)).join('');
    let tsd3 = tsd; if (plantedNoTsdFrom) tsd3 = i < 12 ? (i === 0 ? tsd.slice(0, 5) + (tsd[5] === 'A' ? 'C' : 'A') : tsd) : rb(6);
    let left = rb(25), right = rb(25);
    if (!plantedNoTsdFrom) { if (i === 3) left = left.slice(0, 22) + '-' + left.slice(22); if (i === 4) right = right.slice(0, 2) + '-'.repeat(23); }
    rows.push({ header: 'copy' + i, seq: left + tsd + b + 'AAAAAAAAAA' + tsd3 + right }); }
  return rows;
}
const P = { minLen: 4, maxLen: 20, maxDiv: 20, flankSize: 30 };
const a = T._findTSD(mk(9, 120, 20, true), 'auto', P);
let planted = 0; for (let i = 1; i <= 11; i++) { const r = a.find(x => x.seqName === 'copy' + i); if (r && r.tsdLen === 6 && r.upTSD === r.downTSD) planted++; else console.log('  planted miss', i, r && r.upTSD + '/' + r.downTSD); }
const b = T._findTSD(mk(5, 100, 20, false), 'auto', P);
const mm = b.filter(r => r.mismatches > 0).map(r => r.seqName + ':' + r.upTSD + '/' + r.downTSD);
console.log(process.env.TSDLIB, process.env.TLOOK, process.env.TPEN, process.env.TBASE, `planted ${planted}/11, flank-check pairs with mismatches: ${mm.length} ${mm.join(' ')}`);
