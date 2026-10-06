'use strict';
// ViewAlign's MACSE v2.07 port (va/macse-align.js), vertebrate mitochondrial code (gc 2).
//  main  : per gene, Sb1 own mtDNA (reliable) + the 54 GenBank betulina mitogenomes (less reliable) -> genes/<g>.main_NT/AA
//  pairs : per gene, every other sequence aligned alone with the reference (reference reliable, the other less
//          reliable), so no group of sequences can outvote the reference (82 near-identical nuclear copies did:
//          MACSE placed their shared ND3 indel as a frameshift in every mitogenome) -> genes/<g>.pairs.jsonl
// usage: node run_pairs.js main|pairs GENE[,GENE...]
const fs = require('fs'), path = require('path');
const MA = require('./va/macse-align.js');
const REF = 'Sb1_own_mtDNA_polished';
const mode = process.argv[2], genes = process.argv[3].split(',');
(async () => {
    for (const g of genes) {
        const recs = MA.parseFasta(fs.readFileSync(`genes/${g}.fna`, 'utf8'));
        const ref = recs.find(r => r.name === REF);
        const t = Date.now();
        if (mode === 'main') {
            const set = recs.filter(r => r.name === REF || r.name.startsWith('Sbet_'));   // betulina only: with the 7 other-species
            // mitogenomes in the set, MACSE shifted ND3 of every betulina row over ~75 codons (profile vs profile)
            // all sequences reliable (MACSE default): with the GenBank rows marked less reliable, MACSE (and the port,
            // identically) shifts ND3 of every row over ~75 codons, and the view would show frameshifts that are not there
            const r = MA.alignSequences(set, { gc: 2 });
            fs.writeFileSync(`genes/${g}.main_NT.fna`, MA.toFasta(r.nt));
            fs.writeFileSync(`genes/${g}.main_AA.faa`, MA.toFasta(r.aa));
            console.log(`${g}\tmain\t${set.length} seqs\t${((Date.now() - t) / 1000).toFixed(1)} s`);
        } else {
            const out = [];
            for (const x of recs) {
                if (x.name === REF) continue;
                const r = MA.alignSequences([ref, x], { gc: 2, lessReliable: [x.name] });
                const by = n => r.nt.find(z => z.name === n).seq, byA = n => r.aa.find(z => z.name === n).seq;
                out.push(JSON.stringify({ name: x.name, ref: by(REF), seq: by(x.name), refAA: byA(REF), aa: byA(x.name) }));
            }
            fs.writeFileSync(`genes/${g}.pairs.jsonl`, out.join('\n') + '\n');
            console.log(`${g}\tpairs\t${out.length}\t${((Date.now() - t) / 1000).toFixed(1)} s`);
        }
    }
})().catch(e => { console.error(e); process.exit(1); });
