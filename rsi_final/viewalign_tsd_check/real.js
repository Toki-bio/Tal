// Real rsi plates (C:/work/Tal/rsi_final/stage8_plates): ViewAlign's TSD finder (auto, defaults) against an
// exhaustive oracle with SINEderella's has_tsd rule: 10-20 bp, <= 10 % mismatches, 5' copy ends <= 4 bp before the
// copy's first element residue, 3' copy starts 0-30 bp after its last element residue (element = REF/proposed columns).
// Null: 5' side of copy i with the 3' side of copy i+1.
const fs = require('fs'), path = require('path');
const libs = (process.env.LIBS || './tsdlib.js,./protoD.js').split(',').map(p => [p, require(p)]);
const DIR = process.env.PLATES || __dirname + '/../stage8_plates';
const readFa = f => { const out = []; let cur = null; for (const l of fs.readFileSync(f, 'utf8').split(/\r?\n/)) { if (l[0] === '>') { cur = { header: l.slice(1).trim(), seq: '' }; out.push(cur); } else if (cur) cur.seq += l.trim(); } return out; };
const isGap = c => c === '-' || c === '.';
function oracle(upB, dnB) {   // upB: residues before the element (closest last), dnB: after its end
  let best = null;
  for (let L = 20; L >= 10 && !best; L--) {
    const mmMax = Math.floor(L * 0.1);
    for (let d5 = 0; d5 <= 4 && !best; d5++) {
      if (upB.length < L + d5) continue;
      const a = upB.slice(upB.length - L - d5, upB.length - d5);
      for (let d3 = 0; d3 <= 30 && d3 + L <= dnB.length; d3++) {
        const b = dnB.slice(d3, d3 + L); let mm = 0;
        for (let i = 0; i < L && mm <= mmMax; i++) if (a[i] !== b[i]) mm++;
        if (mm <= mmMax) { best = { L, d5, d3, a, b }; break; }
      }
    }
  }
  return best;
}
const tot = {};
for (const f of fs.readdirSync(DIR).filter(x => x.endsWith('.aln.fa')).sort()) {
  const rows = readFa(path.join(DIR, f));
  const ref = rows.find(r => /^REF_/.test(r.header)), prop = rows.find(r => /^proposed/.test(r.header));
  const copies = rows.filter(r => r !== ref && r !== prop);
  if (!ref || copies.length < 3) continue;
  const span = s => { let a = -1, b = -1; for (let c = 0; c < s.length; c++) if (!isGap(s[c])) { if (a < 0) a = c; b = c; } return [a, b]; };
  const [e0] = span(ref.seq), [, e1] = span(prop ? prop.seq : ref.seq);
  // per copy: residues before the element columns, after them
  const parts = copies.map(r => { const up = [], dn = []; let inEl = 0;
    for (let c = 0; c < r.seq.length; c++) { const ch = r.seq[c].toUpperCase(); if (isGap(ch)) continue; if (c < e0) up.push(ch); else if (c > e1) dn.push(ch); else inEl++; }
    return { up: up.join(''), dn: dn.join(''), inEl }; });
  const orc = parts.map(p => p.inEl > 20 ? oracle(p.up, p.dn) : null);
  const nul = parts.map((p, i) => oracle(p.up, parts[(i + 1) % parts.length].dn));
  const line = [f.replace('.plate.aln.fa', '').padEnd(16), `${copies.length} copies`.padEnd(11), `oracle ${orc.filter(Boolean).length}`, `(null ${nul.filter(Boolean).length})`];
  // control: each copy's 3' side (columns past the element) from the next copy: no TSD can be real there
  const ctrl = rows.map(r => r === ref || r === prop ? r : null);
  copies.forEach((r, i) => { const nx = copies[(i + 1) % copies.length]; const idx = rows.indexOf(r); ctrl[idx] = { header: r.header, seq: r.seq.slice(0, e1 + 1) + nx.seq.slice(e1 + 1) }; });
  for (const [name, T] of libs) {
    const resC = T._findTSD(ctrl, 'auto', { minLen: 4, maxLen: 20, maxDiv: 20, flankSize: 30 });
    const fp10 = resC.filter(r => r.upTSD.length >= 10).length, fp8 = resC.filter(r => r.upTSD.length >= 8).length;
    const res = T._findTSD(rows, 'auto', { minLen: 4, maxLen: 20, maxDiv: 20, flankSize: 30 });
    const by = new Map(res.map(r => [rows[r.seqIndex], r]));
    let agree = 0, missed = 0, other = 0, extra = 0;
    copies.forEach((r, i) => {
      const o = orc[i], v = by.get(r);
      const vseq = v ? v.upTSD : '';
      if (o) { if (!v) missed++; else if (v.upTSD.length >= 6 && (o.a.includes(v.upTSD.slice(0, 6)) || v.upTSD.includes(o.a.slice(0, 6)))) agree++; else other++; }
      else if (v && v.upTSD.length >= 10) extra++;
    });
    tot[name] = tot[name] || { agree: 0, missed: 0, other: 0, extra: 0, o: 0, fp10: 0, fp8: 0, n: 0 };
    tot[name].fp10 += fp10; tot[name].fp8 += fp8; tot[name].n += copies.length;
    Object.assign(tot[name], { agree: tot[name].agree + agree, missed: tot[name].missed + missed, other: tot[name].other + other, extra: tot[name].extra + extra, o: tot[name].o + orc.filter(Boolean).length });
    line.push(`| ${path.basename(name, '.js')}: found ${agree}, missed ${missed}, other ${other}, >=10bp w/o oracle ${extra}  [${res.boundaries.label}]`);
  }
  console.log(line.join(' '));
}
for (const [n, t] of Object.entries(tot)) console.log(`TOTAL ${n}: oracle TSDs ${t.o}; found ${t.agree} (${(100 * t.agree / t.o).toFixed(0)}%), missed ${t.missed}, a different one ${t.other}; >=10 bp calls the oracle has not ${t.extra} | CONTROL (swapped 3' sides, ${t.n} copies): calls >=8 bp ${t.fp8}, >=10 bp ${t.fp10}`);
