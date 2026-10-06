#!/usr/bin/env python3
# Cut the 13 mitochondrial protein-coding genes out of every sequence in seqs/all.fa.
# Query = the gene in Sb1 own mtDNA (unit2; coordinates units/genes_unit2.tsv, lifted from NC_069019).
# blastn per gene; per subject the best HSP plus collinear HSPs on the same strand within the gene's span;
# the extracted region is extended by a query end that did not align only when that end is <= 60 bp (a divergent
# gene end); a longer unaligned end means the copy stops there (NUMT fragment), and its flank is not taken.

import subprocess, os, sys, collections
os.chdir(os.path.expanduser('~/Sicista2026/mito/pseudo_check'))
os.makedirs('genes', exist_ok=True)
def readfa(p):
    d = collections.OrderedDict(); n = None
    for l in open(p):
        l = l.strip()
        if l.startswith('>'): n = l[1:].split()[0]; d[n] = []
        elif n: d[n].append(l)
    return collections.OrderedDict((k, ''.join(v).upper()) for k, v in d.items())
comp = str.maketrans('ACGTRYKMSWBDHVN', 'TGCAYRMKSWVHDBN')
rc = lambda s: s.translate(comp)[::-1]
ext = lambda n: n if n <= 60 else 0
seqs = readfa('seqs/all.fa'); u2 = seqs['Sb1_own_mtDNA_polished']
PCG = ['ND1','ND2','COX1','COX2','ATP8','ATP6','COX3','ND3','ND4L','ND4','ND5','ND6','CYTB']
coords = {}
for l in open('genes_polished.tsv'):
    f = l.rstrip('\n').split('\t')
    if f[0] in PCG: coords[f[0]] = (int(f[1]), int(f[2]))
subprocess.run('makeblastdb -in seqs/all.fa -dbtype nucl -out seqs/alldb > /dev/null', shell=True, check=True)
summ = open('genes/extraction.tsv', 'w')
summ.write('gene\tseq\tcoverage\tbest_pident\tstrand\tstart\tend\tlength\tn_hsp\n')
for g in PCG:
    s, e = coords[g]; q = u2[s-1:e]
    if g == 'ND6': q = rc(q)
    L = len(q)
    open(f'genes/{g}.query.fa', 'w').write(f'>{g}\n{q}\n')
    out = subprocess.run(f'blastn -task blastn -query genes/{g}.query.fa -db seqs/alldb -evalue 1e-5 -dust no '
                         f'-max_target_seqs 5000 -max_hsps 50 -outfmt "6 sseqid pident qstart qend sstart send bitscore"',
                         shell=True, check=True, capture_output=True, text=True).stdout
    hs = collections.defaultdict(list)
    for l in out.splitlines():
        sid, pid, qs, qe, ss, se, bs = l.split('\t')
        hs[sid].append((float(bs), float(pid), int(qs), int(qe), int(ss), int(se)))
    rec = []
    for sid, H in hs.items():
        H.sort(reverse=True); best = H[0]; plus = best[4] <= best[5]
        bs0 = min(best[4], best[5])
        chain = [best]
        for h in H[1:]:
            if (h[4] <= h[5]) != plus: continue
            if abs(min(h[4], h[5]) - bs0) > L + 300: continue
            # collinear with the best HSP: subject offset agrees with the query offset within 300 bp
            dq = h[2] - best[2]; ds = (min(h[4], h[5]) - bs0) if plus else (bs0 - min(h[4], h[5]))
            if abs(dq - ds) > 300: continue
            if any(not (h[3] < c[2] or h[2] > c[3]) and min(h[3], c[3]) - max(h[2], c[2]) > 0.5 * (h[3] - h[2]) for c in chain): continue
            chain.append(h)
        qmin = min(c[2] for c in chain); qmax = max(c[3] for c in chain)
        cov = sum(1 for p in range(1, L+1) if any(c[2] <= p <= c[3] for c in chain)) / L
        S = seqs[sid]
        if plus:
            a = min(min(c[4], c[5]) for c in chain) - ext(qmin - 1); b = max(max(c[4], c[5]) for c in chain) + ext(L - qmax)
            a = max(a, 1); b = min(b, len(S)); sub = S[a-1:b]
        else:
            a = min(min(c[4], c[5]) for c in chain) - ext(L - qmax); b = max(max(c[4], c[5]) for c in chain) + ext(qmin - 1)
            a = max(a, 1); b = min(b, len(S)); sub = rc(S[a-1:b])
        rec.append((sid, cov, best[1], '+' if plus else '-', a, b, sub, len(chain)))
    order = {n: i for i, n in enumerate(seqs)}
    rec.sort(key=lambda r: order[r[0]])
    with open(f'genes/{g}.fna', 'w') as fo:
        for r in rec:
            summ.write(f'{g}\t{r[0]}\t{r[1]:.3f}\t{r[2]:.1f}\t{r[3]}\t{r[4]}\t{r[5]}\t{len(r[6])}\t{r[7]}\n')
            if r[1] >= 0.5: fo.write(f'>{r[0]}\n{r[6]}\n')
    print(g, L, 'hits', len(rec), 'kept', sum(1 for r in rec if r[1] >= 0.5), 'missing', len(seqs) - len(rec))
print('GENES_DONE')
