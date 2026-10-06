'use strict';
// ViewAlign check on real data (Sicista mitochondrial genes, 2026-10-06): load the result alignments with their gene BED,
// turn on codon analysis (by gene), and compare per row what the viewer computed with an independent count in Node:
//  stops      : stop codons per row (internal + terminal, incl. polyA-completed T/TA), translating each gene in the
//               row's own phase (gaps skipped), vertebrate mito code
//  syn / non  : bases marked synonymous / non-synonymous against row 1 (ViewAlign's rule: each base that differs from
//               the reference codon starting at the same column, classed by that single change)
//  frameshift : rows the viewer marks with any frameshift vs rows the pipeline (pairwise MACSE) says have one
// Also saves screenshots for the manual-inspection guide.
const fs = require('fs'), path = require('path');
const REPO = 'C:/work/MSA-viewer';
const { start } = require(path.join(REPO, 'tests/lib/static-server'));
const { launch, loadFasta } = require(path.join(REPO, 'tests/lib/browser'));
const DIR = 'C:/work/hylomys_ont/sicista_mito_pseudo';
const OUT = process.argv[2] || DIR + '/viewalign_check';
fs.mkdirSync(OUT, { recursive: true });

const B4 = 'TCAG', AAS = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
const CODE = {}; let k = 0;
for (const a of B4) for (const b of B4) for (const c of B4) CODE[a + b + c] = AAS[k++];
Object.assign(CODE, { TGA: 'W', ATA: 'M', AGA: '*', AGG: '*' });
const parseFa = t => { const r = []; for (const blk of t.split('>').slice(1)) { const i = blk.indexOf('\n'); r.push({ name: blk.slice(0, i).trim(), seq: blk.slice(i + 1).replace(/\s/g, '') }); } return r; };

function expected(rows, bed) {
    // BED positions are positions in row 1's own sequence; map them to alignment columns through its gaps
    const colOf = []; [...rows[0].seq].forEach((ch, c) => { if (ch !== '-' && ch !== '.') colOf.push(c); });
    const genes = bed.trim().split('\n').map(l => { const f = l.split('\t'); return { s: colOf[+f[1]], e: colOf[+f[2] - 1] + 1 }; });
    return rows.map((r, ri) => {
        let stops = 0, syn = 0, non = 0;
        for (const g of genes) {
            const L = g.e - g.s;
            const sub = rows.map(x => x.seq.slice(g.s, g.e));
            // codons of a row in its own phase: [{col, codon}]
            const cods = s => {
                const out = []; let buf = '', cols = [], seen = false;
                for (let c = 0; c < L; c++) {
                    const ch = s[c];
                    if (ch === '-' || ch === '.') continue;
                    if (!seen) { seen = true; const lead = c % 3; if (lead) { buf = '?'.repeat(lead); } }
                    buf += ch.toUpperCase().replace('U', 'T'); cols.push(c);
                    if (buf.length === 3) { if (!buf.includes('?')) out.push({ col: cols[0], codon: buf }); buf = ''; cols = []; }
                }
                return { out, tail: buf };
            };
            const mine = cods(sub[ri]);
            for (const x of mine.out) if (CODE[x.codon] === '*') stops++;
            if (mine.tail === 'T' || mine.tail === 'TA') stops++;
            if (ri === 0) continue;
            const ref = new Map(cods(sub[0]).out.map(x => [x.col, x.codon]));
            for (const x of mine.out) {
                const rc = ref.get(x.col);
                if (!rc || rc === x.codon || !/^[ACGT]{3}$/.test(rc) || !/^[ACGT]{3}$/.test(x.codon)) continue;
                const aa = CODE[x.codon], raa = CODE[rc];
                for (let p = 0; p < 3; p++) {
                    if (x.codon[p] === rc[p]) continue;
                    if (aa === raa) { syn++; continue; }
                    const m = rc.slice(0, p) + x.codon[p] + rc.slice(p + 1);
                    if ((CODE[m] || 'X') !== raa) non++; else syn++;
                }
            }
        }
        return { stops, syn, non };
    });
}

(async () => {
    const { server, baseUrl } = await start();
    const browser = await launch();
    const page = await browser.newPage({ viewport: { width: 1600, height: 1000 } });
    const report = [];
    for (const set of ['mitogenomes', 'all', 'numtB_family']) {
        const fasta = fs.readFileSync(`${DIR}/viewer_${set}_13genes.fasta`, 'utf8');
        const bed = fs.readFileSync(`${DIR}/viewer_${set}_13genes.bed`, 'utf8');
        const rows = parseFa(fasta);
        await page.goto(baseUrl + '/index.html', { waitUntil: 'networkidle' });
        await loadFasta(page, fasta);
        const t0 = Date.now();
        const v = await page.evaluate((bed) => {
            setAnnotation(bed, 'genes.bed');
            document.getElementById('codonFrame').value = 'auto';
            const cb = document.getElementById('codonAnalysis'); cb.checked = true; cb.dispatchEvent(new Event('change'));
            return new Promise(res => setTimeout(() => {
                const cd = state._codonData;
                res({
                    code: document.getElementById('codonCode').value, byGene: !!(cd && cd.byGene), genes: cd && cd.genes,
                    rows: state.seqs.map((s, i) => ({
                        name: s.fullHeader || s.header,
                        stops: cd.aaSeq[i].filter(e => e.aa === '*').length,
                        fs: cd.frameShifts[i].length,
                        syn: cd.synNonSyn[i].filter(x => x === 'syn').length,
                        non: cd.synNonSyn[i].filter(x => x === 'nonsyn').length,
                    })),
                });
            }, 1500));
        }, bed);
        const exp = expected(rows, bed);
        const pipe = Object.fromEntries(fs.readFileSync(`${DIR}/per_sequence.tsv`, 'utf8').trim().split('\n').slice(1).map(l => { const f = l.split('\t'); return [f[0], { fs: +f[4], stops: +f[5] }]; }));
        let bad = [];
        const fsTab = { both: 0, viewerOnly: 0, pipelineOnly: 0, neither: 0 };
        v.rows.forEach((r, i) => {
            const e = exp[i];
            if (r.name !== rows[i].name) bad.push(`row ${i} name ${r.name} vs ${rows[i].name}`);
            if (r.stops !== e.stops || r.syn !== e.syn || r.non !== e.non) bad.push(`${r.name.split('|')[0]}: viewer stops/syn/non ${r.stops}/${r.syn}/${r.non}, expected ${e.stops}/${e.syn}/${e.non}`);
            const p = pipe[r.name.split('|')[0]];
            const pv = p && p.fs > 0, vv = r.fs > 0;
            fsTab[pv && vv ? 'both' : vv ? 'viewerOnly' : pv ? 'pipelineOnly' : 'neither']++;
        });
        report.push(`## ${set}: ${rows.length} rows, ${rows[0].seq.length} columns; code ${v.code}, by gene ${v.byGene}, genes ${v.genes}; ${((Date.now() - t0) / 1000).toFixed(1)} s`);
        report.push(`stops/syn/non per row, viewer vs independent count: ${bad.length ? bad.length + ' rows differ' : 'all ' + rows.length + ' rows identical'}`);
        bad.slice(0, 15).forEach(b => report.push('  ' + b));
        report.push(`frameshift rows (viewer any mark vs pipeline persistent shift): ${JSON.stringify(fsTab)}`);
        fs.writeFileSync(`${OUT}/viewer_rows_${set}.tsv`, 'row\tstops\tframeshift_marks\tsyn\tnonsyn\texp_stops\texp_syn\texp_non\n' +
            v.rows.map((r, i) => `${r.name}\t${r.stops}\t${r.fs}\t${r.syn}\t${r.non}\t${exp[i].stops}\t${exp[i].syn}\t${exp[i].non}`).join('\n') + '\n');
        await page.screenshot({ path: `${OUT}/screen_${set}.png` });
    }
    fs.writeFileSync(`${OUT}/REPORT.txt`, report.join('\n') + '\n');
    console.log(report.join('\n'));
    await browser.close(); server.close();
})().catch(e => { console.error(e); process.exit(1); });
