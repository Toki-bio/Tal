// copies with no TSD (random 3' side): how many get a TSD call anyway, by length
const T = require(process.env.TSDLIB);
const { build, setSeed } = require('./sim.js');
let n = 0, calls = 0; const byLen = {};
for (let k = 0; k < 20; k++) { setSeed(1000 + k * 7919); const rows = build({ noTsd: true });
  const res = T._findTSD(rows, 'auto', { minLen: 4, maxLen: 20, maxDiv: 20, flankSize: 30 });
  n += rows.length; res.forEach(r => { calls++; const L = r.upPositions.length; const b = L >= 8 ? '>=8' : L >= 6 ? '6-7' : '4-5'; byLen[b] = (byLen[b] || 0) + 1; }); }
console.log(`${process.env.TSDLIB}: TSD called on ${(100 * calls / n).toFixed(0)}% of copies without one`, JSON.stringify(byLen));
