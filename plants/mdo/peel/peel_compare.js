// Profiles the reconstructed groups on the 496-seed plate (the measures of SINE_discriminator RSI_PEEL_LOG.md)
// and scores ViewAlign's existing peel and tree cut (MSA-viewer peel.js, kmer-tree.js) against them.
//   node peel_compare.js [path to the ViewAlign checkout]
// Reads ../alignments/mdo_consensuses.aln.fa and seeds_vs_peeled.tsv; writes peel_compare.txt.
// Group membership is a reconstruction: best consensus with identity >= 0.6 over >= 0.5 of the seed.
const fs = require('fs'), path = require('path');
const VA = process.argv[2] || 'C:/work/MSA-viewer-new';
const KT = require(path.join(VA, 'kmer-tree.js')), Peel = require(path.join(VA, 'peel.js'));
const { ari } = require(path.join(VA, 'tests/kmer/metrics.js'));
const H = __dirname;
const seqs = fs.readFileSync(path.join(H, '../alignments/mdo_consensuses.aln.fa'), 'utf8').split('>').slice(1)
  .map(b => { const l = b.split('\n'); return { header: l[0].trim(), seq: l.slice(1).join('').trim().toUpperCase() }; });
const tsv = fs.readFileSync(path.join(H, 'seeds_vs_peeled.tsv'), 'utf8').trim().split('\n').slice(1).map(l => l.split('\t'));
const grp = new Map(tsv.map(f => [f[0], (+f[4] >= 0.6 && +f[5] >= 0.5) ? f[2] : null]));
const truth = seqs.map(s => grp.get(s.header));
const n = seqs.length, L = seqs[0].seq.length;
const out = [];
const say = s => { out.push(s); console.log(s); };

function ident(a, b) {
  let same = 0, both = 0;
  for (let i = 0; i < L; i++) if (a[i] !== '-' && b[i] !== '-') { both++; if (a[i] === b[i]) same++; }
  return both >= 50 ? same / both : null;
}
const sim = seqs.map(a => seqs.map(b => ident(a.seq, b.seq)));
const mean = v => v.length ? v.reduce((x, y) => x + y, 0) / v.length : NaN;

say(`plate ${n} sequences, ${L} columns; in a reconstructed group ${truth.filter(Boolean).length}, outside ${truth.filter(t => !t).length}`);
say('group\tsize\twithin\tvs_rest\tdiag_cols\tof_which_gap_state\tplate_rows_first-last\tlongest_adjacent_run');
const names = [...new Set(truth.filter(Boolean))].sort((a, b) => parseInt(a.split('_')[1], 10) - parseInt(b.split('_')[1], 10));
for (const g of names) {
  const m = [], r = [];
  truth.forEach((t, i) => (t === g ? m : r).push(i));
  const within = [], vs = [];
  for (let i = 0; i < m.length; i++) {
    for (let j = i + 1; j < m.length; j++) if (sim[m[i]][m[j]] !== null) within.push(sim[m[i]][m[j]]);
    for (const k of r) if (sim[m[i]][k] !== null) vs.push(sim[m[i]][k]);
  }
  let diag = 0, gaps = 0;
  for (let c = 0; c < L; c++) {
    const count = {};
    for (const i of m) { const ch = seqs[i].seq[c]; count[ch] = (count[ch] || 0) + 1; }
    const [ch, k] = Object.entries(count).sort((a, b) => b[1] - a[1])[0];
    if (k / m.length >= 0.9 && r.filter(i => seqs[i].seq[c] === ch).length / r.length <= 0.5) { diag++; if (ch === '-') gaps++; }
  }
  let run = 0, best = 0;
  truth.forEach(t => { run = t === g ? run + 1 : 0; best = Math.max(best, run); });
  say([g, m.length, mean(within).toFixed(3), mean(vs).toFixed(3), diag, gaps, `${m[0] + 1}-${m[m.length - 1] + 1}`, best].join('\t'));
}

const idx = truth.map((t, i) => (t ? i : -1)).filter(i => i >= 0);
const lab = (groups, un) => { const o = new Array(n).fill(null); groups.forEach((g, gi) => g.forEach(i => { o[i] = 'g' + gi; })); un.forEach(i => { o[i] = 'u' + i; }); return idx.map(i => o[i]); };
const tr = idx.map(i => truth[i]);
const tree = KT.guideTree(seqs, 6, { metric: 'pdist' });
say('\nmethod (ViewAlign code)\tgroups\tARI_against_reconstructed_groups\tunassigned');
for (const [label, k] of [['tree cut, automatic', 'auto'], [`tree cut, ${names.length} groups given`, names.length]]) {
  const c = KT.cutTree(tree, k, 3);
  say([label, c.groups.length, ari(tr, lab(c.groups, c.unassigned)).toFixed(3), c.unassigned.length].join('\t'));
}
for (const outliers of [false, true]) {
  const p = Peel.peel(seqs, { metric: 'pdist', minSize: 3, minGap: 0.04, outlierZ: 10, outliers });
  say([`peel loop, outliers ${outliers ? 'on' : 'off'}`, p.groups.length, ari(tr, lab(p.groups, p.unassigned)).toFixed(3), p.unassigned.length].join('\t'));
}
fs.writeFileSync(path.join(H, 'peel_compare.txt'), out.join('\n') + '\n');
