// Planted-TSD alignments run through ViewAlign's own _findTSD (auto mode, default parameters).
// Each copy: 5' flank | TSD | SINE body (mutated) | poly-A tail | TSD (0-1 mismatch) | 3' flank.
// Columns: flank block, body block, tail block (gaps padded), flank block - as an MSA of copies with flanks looks.
const T = require(process.env.TSDLIB || './tsdlib.js');
let seed = 1;
const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff; };
const ri = (a, b) => a + Math.floor(rnd() * (b - a + 1));
const base = () => 'ACGT'[ri(0, 3)];
const randSeq = n => Array.from({ length: n }, base).join('');
const mutate = (s, p) => s.split('').map(c => rnd() < p ? 'ACGT'.replace(c, '')[ri(0, 2)] : c).join('');

const DEF = { n: 20, body: 180, div: 0.12, tail: [12, 12], tsd: [8, 16], tsdMis: 0.3, flank: 50, trunc5: 0, trunc3: 0, tsdIndel: 0, firstGappy: false, tailGapsAt: 'end', aTsd: 0 };
function build(o) {
  o = Object.assign({}, DEF, o);
  const cons = randSeq(o.body);
  const rows = [];
  const maxTail = o.tail[1];
  for (let i = 0; i < o.n; i++) {
    let tsd = randSeq(ri(o.tsd[0], o.tsd[1]));
    if (rnd() < o.aTsd) tsd = 'AAAA' + tsd.slice(4);
    let tsd2 = rnd() < o.tsdMis ? mutate(tsd, 1 / tsd.length) : tsd;
    if (o.noTsd) tsd2 = randSeq(tsd.length);
    if (rnd() < o.tsdIndel) { const p = ri(2, tsd2.length - 3); tsd2 = tsd2.slice(0, p) + tsd2.slice(p + 1); }
    let body = mutate(cons, o.div);
    let lead = '';
    if (rnd() < o.trunc5) { const cut = ri(20, 80); lead = '-'.repeat(cut); body = body.slice(cut); }
    const tl = ri(o.tail[0], o.tail[1]);
    let tail = 'A'.repeat(tl);
    tail = o.tailGapsAt === 'end' ? tail + '-'.repeat(maxTail - tl) : '-'.repeat(maxTail - tl) + tail;
    const f5 = randSeq(o.flank), f3 = randSeq(o.flank);
    // 5' flank block: flank + TSD right-justified against the body (fixed width)
    const left = (f5 + tsd).padStart(o.flank + o.tsd[1], '-');
    const right = (tsd2 + f3).padEnd(o.flank + o.tsd[1], '-');
    const seq = left + lead + body + tail + right;
    rows.push({ header: 'copy' + (i + 1), seq, tsd, tsd2, tsdUpEndCol: left.length - 1, tsdDownStartCol: left.length + lead.length + body.length + tail.length });
  }
  if (o.firstGappy) { const r = rows[0]; r.seq = r.seq.slice(0, o.flank + o.tsd[1] + 30) + '-'.repeat(25) + r.seq.slice(o.flank + o.tsd[1] + 55); }
  return rows;
}
function score(rows, params = {}) {
  const p = Object.assign({ minLen: 4, maxLen: 20, maxDiv: 20, flankSize: 30 }, params);
  const res = T._findTSD(rows, 'auto', p);
  let hit = 0, wrong = 0, miss = 0, skipped = 0;
  const found = new Map(res.map(r => [r.seqIndex, r]));
  rows.forEach((row, i) => {
    const r = found.get(i);
    const inMiss = res.misses.some(m => m.seqIndex === i);
    if (!r) { if (inMiss) miss++; else skipped++; return; }
    // right if the up copy ends at the planted up TSD end and the down copy starts at the planted down start (+-1)
    // right: covers >= 80% of the planted copy on both sides and is at most 4 bases longer (a longer exact
    // duplication is the same call by sequence)
    const L = row.tsd.length, upSet = new Set(Array.from({ length: L }, (_, j) => row.tsdUpEndCol - j));
    const downSet = new Set(Array.from({ length: row.tsd2.length }, (_, j) => row.tsdDownStartCol + j));
    const cu = r.upPositions.filter(c => upSet.has(c)).length, cd = r.downPositions.filter(c => downSet.has(c)).length;
    if (cu >= 0.8 * L && cd >= 0.8 * row.tsd2.length && r.upPositions.length <= L + 4) hit++; else wrong++;
  });
  return { hit, wrong, miss, skipped, boundary: res.boundaries.label };
}
const scen = {
  'baseline: equal tails (12 A)': {},
  'poly-A tail 8-30, gaps after the A run': { tail: [8, 30] },
  'poly-A tail 8-30, gaps before the A run': { tail: [8, 30], tailGapsAt: 'start' },
  '30% copies 5\'-truncated': { trunc5: 0.3 },
  'TSD with a 1-bp indel in 30%': { tsdIndel: 0.3 },
  'old copies, 20% body divergence': { div: 0.2 },
  'TSD starting with AAAA (A-rich target) in 50%': { aTsd: 0.5 },
  'first row a real copy with a 25-bp gap': { firstGappy: true },
  'only 5 copies': { n: 5 },
  'TSD 15-25 bp (max length 20)': { tsd: [15, 25] },
};
const reps = 20;
if (require.main === module) for (const [name, o] of Object.entries(scen)) {
  const tot = { hit: 0, wrong: 0, miss: 0, skipped: 0 };
  let lab = '';
  for (let k = 0; k < reps; k++) { seed = 1000 + k * 7919; const s = score(build(o)); for (const x in tot) tot[x] += s[x]; lab = s.boundary; }
  const n = reps * (o.n || DEF.n);
  console.log(`${name.padEnd(48)} right ${(100 * tot.hit / n).toFixed(0).padStart(3)}%  wrong ${(100 * tot.wrong / n).toFixed(0).padStart(3)}%  missed ${(100 * tot.miss / n).toFixed(0).padStart(3)}%  skipped ${(100 * tot.skipped / n).toFixed(0).padStart(3)}%   (${lab})`);
}
module.exports = { build, score, setSeed: v => { seed = v; } };
