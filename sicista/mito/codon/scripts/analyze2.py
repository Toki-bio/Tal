#!/usr/bin/env python3
# Pseudogene signals in Sicista mitochondrial protein-coding genes (2026-10-06), from pairwise MACSE (ViewAlign port)
# alignments of every sequence with the read-polished own mtDNA (reliable; gc 2). Everything is measured in the
# reference's codon coordinates, so no group of sequences can shift another's frame.
#  frameshifts : a run where the sequence's frame stays shifted against the reference for > 30 columns (10 codons) and
#                starts before the last 30 columns of the gene (an insertion + deletion a few columns apart is MACSE
#                aligning a divergent stretch; gene-end ambiguity is not a frameshift of the gene)
#  stops       : internal '*' in MACSE's translation of the sequence (terminal stop excluded)
#  pN/pS       : each member against its group's consensus codon, Nei-Gojobori counting, vertebrate mito code;
#                functional mtDNA within a species ~0.1, a neutral pseudogene ~1
#  lineage     : betulina GenBank records closer to own mtDNA (A) or to the nuclear B-family consensus (B)
#  NUMT signal : NUMT-private alleles carried (positions where the nuclear B family differs from both GenBank
#                lineage consensus sequences) and N placement at positions where the nuclear family differs
import os, re, json, math, collections, itertools
os.chdir(os.path.expanduser('~/Sicista2026/mito/pseudo_check'))
PCG = ['ND1', 'ND2', 'COX1', 'COX2', 'ATP8', 'ATP6', 'COX3', 'ND3', 'ND4L', 'ND4', 'ND5', 'ND6', 'CYTB']
REF = 'Sb1_own_mtDNA_polished'
B4 = 'TCAG'
AAS = 'FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODE = {a + b + c: AAS[16 * i + 4 * j + k] for i, a in enumerate(B4) for j, b in enumerate(B4) for k, c in enumerate(B4)}
CODE.update({'TGA': 'W', 'ATA': 'M', 'AGA': '*', 'AGG': '*'})
GROUP_ORDER = ['REF_own_mtDNA', 'GB_betulina_A', 'GB_betulina_B', 'GB_other_species',
               'NUMT_A_young', 'NUMT_B_family', 'NUMT_old', 'NUMT_trizona']
ACGT3 = re.compile('[ACGT]{3}')


def group0(n):
    if n == REF: return 'REF_own_mtDNA'
    for pre, gname in [('Sb1_NUMTB', 'NUMT_B_family'), ('Sb1_NUMTA', 'NUMT_A_young'), ('Sb1_NUMTold', 'NUMT_old'),
                       ('Stri_NUMT', 'NUMT_trizona'), ('Sbet_', 'GB_betulina')]:
        if n.startswith(pre): return gname
    return 'GB_other_species'


SITES = {}
for c, a in CODE.items():
    if a != '*':
        s = sum(sum(1 for x in 'ACGT' if x != c[p] and CODE[c[:p] + x + c[p + 1:]] == a) / 3 for p in range(3))
        SITES[c] = (s, 3 - s)


def ng_diff(c1, c2):
    d = [p for p in range(3) if c1[p] != c2[p]]
    if not d: return 0.0, 0.0
    sd = nd = 0.0; npath = 0
    for order in itertools.permutations(d):
        cur = c1; s = n = 0; bad = False
        for p in order:
            nxt = cur[:p] + c2[p] + cur[p + 1:]
            if CODE[nxt] == '*' and nxt != c2: bad = True; break
            if CODE[cur] == CODE[nxt]: s += 1
            else: n += 1
            cur = nxt
        if bad: continue
        sd += s; nd += n; npath += 1
    return (sd / npath, nd / npath) if npath else (0.0, float(len(d)))


refseq = {g: [l.strip() for l in open(f'genes/{g}.query.fa')][1] for g in PCG}
names = [l[1:].split()[0] for l in open('seqs/all.fa') if l.startswith('>')]
proj, ins, fsn, stops, ncount = {}, {}, {}, {}, {}
fs_where, stop_where = collections.defaultdict(list), collections.defaultdict(list)
for g in PCG:
    Lg = len(refseq[g])
    proj[g] = {REF: list(refseq[g])}
    ins[g] = {REF: [''] * (Lg + 1)}
    for line in open(f'genes/{g}.pairs.jsonl'):
        r = json.loads(line); n = r['name']; R = r['ref']; S = r['seq']
        p = ['-'] * Lg; insl = [''] * (Lg + 1); rp = 0
        idx = [i for i, ch in enumerate(S) if ch not in '-!']
        a0, a1 = (idx[0], idx[-1]) if idx else (1, 0)
        off = 0; start = None; nfs = 0
        for i, (rc, sc) in enumerate(zip(R, S)):
            if rc not in '-!':
                p[rp] = sc if sc not in '!' else '-'
                rp += 1
            elif sc not in '-!':
                insl[rp] += sc
            if a0 <= i <= a1:
                d = 1 if (rc == '!' and sc not in '-!') else -1 if (sc == '!' and rc not in '-!') else 0
                if d:
                    was = off % 3; off += d
                    if was == 0 and off % 3: start = i
                    elif was and off % 3 == 0:
                        if i - start > 30 and start < a1 - 30: nfs += 1; fs_where[n].append(f'{g}@nt{len(R[:start].replace("-", "").replace("!", "")) + 1}')
                        start = None
        if start is not None and off % 3 and start < a1 - 30: nfs += 1; fs_where[n].append(f'{g}@nt{len(R[:start].replace("-", "").replace("!", "")) + 1}-end')
        # stops: premature stop codons at homologous codon positions, read in the reference frame from the projected
        # bases (not from MACSE's translation: with a less-reliable sequence MACSE writes '!' columns into BOTH rows of
        # some pairs (ND2, ND3) and its translation of that stretch is out of frame for both)
        rs = refseq[g]
        st = [k for k in range(Lg // 3) if ACGT3.fullmatch(''.join(p[3 * k:3 * k + 3]).upper())
              and CODE[''.join(p[3 * k:3 * k + 3]).upper()] == '*' and CODE.get(rs[3 * k:3 * k + 3]) != '*']
        if st: stop_where[n].append(f'{g}:codon' + '/'.join(str(k + 1) for k in st))
        proj[g][n] = p; ins[g][n] = insl
        fsn[(g, n)] = nfs; stops[(g, n)] = len(st); ncount[(g, n)] = S.upper().count('N')


def cod(g, n):
    p = proj[g][n]
    return [''.join(p[i:i + 3]).upper() for i in range(0, len(p) - len(p) % 3, 3)]


present = {n: [g for g in PCG if n in proj[g]] for n in names}
per = collections.OrderedDict()
for n in names:
    gs = present[n]
    cov = sum(sum(1 for c in cod(g, n) if ACGT3.fullmatch(c)) for g in gs) / sum(len(refseq[g]) // 3 for g in PCG)
    per[n] = dict(group=group0(n), genes=len(gs), cov=cov, fs=sum(fsn.get((g, n), 0) for g in gs),
                  stops=sum(stops.get((g, n), 0) for g in gs), N=sum(ncount.get((g, n), 0) for g in gs))


def consensus(members, g):
    rows = [cod(g, m) for m in members if m in proj[g]]
    out = []
    for col in zip(*rows):
        cc = collections.Counter(c for c in col if ACGT3.fullmatch(c))
        out.append(cc.most_common(1)[0][0] if cc else '---')
    return out


groups = collections.defaultdict(list)
for n, r in per.items(): groups[r['group']].append(n)
numtb = {g: consensus(groups['NUMT_B_family'], g) for g in PCG}
refc = {g: cod(g, REF) for g in PCG}


def pdist(n, other):
    d = t = 0
    for g in present[n]:
        for a_, b_ in zip(cod(g, n), other[g]):
            for p in range(3):
                if a_[p] in 'ACGT' and b_[p] in 'ACGT': t += 1; d += a_[p] != b_[p]
    return d / t if t else float('nan')


for n in groups['GB_betulina']:
    per[n]['d_own'] = pdist(n, refc); per[n]['d_numtb'] = pdist(n, numtb)
    per[n]['group'] = 'GB_betulina_A' if per[n]['d_own'] < per[n]['d_numtb'] else 'GB_betulina_B'
groups = collections.defaultdict(list)
for n, r in per.items(): groups[r['group']].append(n)
consA = {g: consensus(groups['GB_betulina_A'], g) for g in PCG}
consB = {g: consensus(groups['GB_betulina_B'], g) for g in PCG}


def pnps_pairs(pairs):
    Sd = Nd = S = N = 0.0; ts = tv = stopch = 0
    for c1, c2 in pairs:
        if not (ACGT3.fullmatch(c1) and ACGT3.fullmatch(c2)): continue
        if CODE[c1] == '*' or CODE[c2] == '*':
            stopch += c1 != c2; continue
        S += (SITES[c1][0] + SITES[c2][0]) / 2; N += (SITES[c1][1] + SITES[c2][1]) / 2
        if c1 != c2:
            sd, nd = ng_diff(c1, c2); Sd += sd; Nd += nd
            for p in range(3):
                if c1[p] != c2[p]:
                    if {c1[p], c2[p]} in ({'A', 'G'}, {'C', 'T'}): ts += 1
                    else: tv += 1
    pS = Sd / S if S else float('nan'); pN = Nd / N if N else float('nan')
    return dict(syn=Sd, non=Nd, pS=pS, pN=pN, ratio=pN / pS if pS else float('nan'), ts=ts, tv=tv, stopch=stopch)


def group_pnps(members):
    pairs = []
    for g in PCG:
        mem = [m for m in members if m in proj[g]]
        if len(mem) < 2: continue
        cons = consensus(mem, g)
        for m in mem: pairs += list(zip(cons, cod(g, m)))
    return pnps_pairs(pairs)


gst = {k: group_pnps(v) for k, v in groups.items() if len(v) >= 2}
btw = []
for a_n, a_c, b_n, b_c in [('GenBank A consensus', consA, 'GenBank B consensus', consB),
                           ('GenBank B consensus', consB, 'nuclear B-family consensus', numtb),
                           ('GenBank A consensus', consA, 'nuclear B-family consensus', numtb),
                           ('own mtDNA', refc, 'GenBank A consensus', consA)]:
    s = pnps_pairs([pr for g in PCG for pr in zip(a_c[g], b_c[g])])
    btw.append(f"{a_n} vs {b_n}: syn {s['syn']:.1f}, nonsyn {s['non']:.1f}, pS {s['pS']:.4f}, pN {s['pN']:.4f}, pN/pS {s['ratio']:.3f}")

# mosaic scan on positions where GenBank A and B consensus differ
diag = [(g, i, p, a_[p], b_[p]) for g in PCG for i, (a_, b_) in enumerate(zip(consA[g], consB[g]))
        for p in range(3) if a_[p] in 'ACGT' and b_[p] in 'ACGT' and a_[p] != b_[p]]
mosaic = {}
for n in groups['GB_betulina_A'] + groups['GB_betulina_B']:
    cc = {g: cod(g, n) for g in present[n]}
    calls = ['A' if cc[g][i][p] == a_ else 'B' if cc[g][i][p] == b_ else '?' for g, i, p, a_, b_ in diag if g in cc]
    wins = []
    for k in range(0, len(calls) - 19, 20):
        c = collections.Counter(calls[k:k + 20]); wins.append('A' if c['A'] > c['B'] else 'B')
    mosaic[n] = (len(calls), calls.count('B') / max(1, len(calls)), ''.join(wins), sum(x != y for x, y in zip(wins, wins[1:])))
# NUMT-private alleles
priv = [(g, i, p, n_[p]) for g in PCG for i, (a_, b_, n_) in enumerate(zip(consA[g], consB[g], numtb[g]))
        for p in range(3) if all(x[p] in 'ACGT' for x in (a_, b_, n_)) and n_[p] != a_[p] and n_[p] != b_[p]]
numt_alleles = {}
for n in groups['GB_betulina_A'] + groups['GB_betulina_B'] + [REF]:
    cc = {g: cod(g, n) for g in present[n]}
    hits = [g for g, i, p, x in priv if g in cc and cc[g][i][p] == x]
    numt_alleles[n] = (len(hits), ' '.join(f'{g}:{c}' for g, c in collections.Counter(hits).items()))


def n_enrich(members, lc, test):
    nN = nNd = sites = dsites = 0
    for n in members:
        for g in present[n]:
            for i, c in enumerate(cod(g, n)):
                for p in range(3):
                    if lc[g][i][p] not in 'ACGT' or c[p] == '-': continue
                    t = test(g, i, p)
                    if t is None: continue
                    sites += 1; dsites += t
                    if c[p] == 'N': nN += 1; nNd += t
    exp = nN * dsites / sites if sites else 0
    pv = (1 - sum(math.exp(-exp) * exp ** k / math.factorial(k) for k in range(nNd))) if nNd else 1.0
    return nN, nNd, exp, max(pv, 0.0)


varB = {}
for g in PCG:
    rows = [cod(g, m) for m in groups['GB_betulina_B'] if m in proj[g]]
    for i in range(len(consB[g])):
        for p in range(3):
            varB[(g, i, p)] = int(len({r[i][p] for r in rows if r[i][p] in 'ACGT'}) > 1)
ntest = [
    ('lineage-B records: N at sites where nuclear family differs from B consensus', groups['GB_betulina_B'], consB,
     lambda g, i, p: None if numtb[g][i][p] not in 'ACGT' else int(numtb[g][i][p] != consB[g][i][p])),
    ('lineage-B records: same, only sites invariant within GenBank B', groups['GB_betulina_B'], consB,
     lambda g, i, p: None if numtb[g][i][p] not in 'ACGT' or varB[(g, i, p)] else int(numtb[g][i][p] != consB[g][i][p])),
    ('lineage-B records, control: N at sites variable within GenBank B', groups['GB_betulina_B'], consB,
     lambda g, i, p: varB[(g, i, p)]),
    ('lineage-A records: N at sites where nuclear family differs from A consensus', groups['GB_betulina_A'], consA,
     lambda g, i, p: None if numtb[g][i][p] not in 'ACGT' else int(numtb[g][i][p] != consA[g][i][p])),
]

os.makedirs('results', exist_ok=True)
f4 = lambda x: f'{x:.4f}' if isinstance(x, float) else ''
with open('results/per_sequence.tsv', 'w') as o:
    o.write('sequence\tgroup\tgenes\tcodon_coverage\tframeshifts\tinternal_stops\tN_in_genes\tframeshift_where(gene@reference_nt)\tstop_where(gene:codon)\t'
            'd_to_own_mtDNA\td_to_NUMTB_cons\tNUMT_private_alleles\tNUMT_alleles_by_gene\tmosaic_sites\tmosaic_fracB\tmosaic_windows\tmosaic_switches\n')
    for n, r in per.items():
        mo = mosaic.get(n, ('', '', '', '')); na = numt_alleles.get(n, ('', ''))
        o.write(f"{n}\t{r['group']}\t{r['genes']}\t{r['cov']:.3f}\t{r['fs']}\t{r['stops']}\t{r['N']}\t{','.join(fs_where[n])}\t"
                f"{','.join(stop_where[n])}\t{f4(r.get('d_own'))}\t{f4(r.get('d_numtb'))}\t{na[0]}\t{na[1]}\t{mo[0]}\t{f4(mo[1])}\t{mo[2]}\t{mo[3]}\n")
with open('results/groups.tsv', 'w') as o:
    o.write('group\tn_seqs\tframeshifts\tseqs_with_frameshift\tinternal_stops\tseqs_with_stop\tsyn_diff\tnonsyn_diff\tpS\tpN\tpN_pS\tts\ttv\tstop_changes_vs_consensus\n')
    for k in GROUP_ORDER:
        v = groups.get(k, [])
        if not v: continue
        s = gst.get(k)
        tail = (f"{s['syn']:.1f}\t{s['non']:.1f}\t{s['pS']:.4f}\t{s['pN']:.4f}\t{s['ratio']:.3f}\t{s['ts']}\t{s['tv']}\t{s['stopch']}" if s else '\t' * 7)
        o.write(f"{k}\t{len(v)}\t{sum(per[n]['fs'] for n in v)}\t{sum(1 for n in v if per[n]['fs'])}\t"
                f"{sum(per[n]['stops'] for n in v)}\t{sum(1 for n in v if per[n]['stops'])}\t{tail}\n")
    o.write('\n# between consensus sequences\n' + ''.join('# ' + x + '\n' for x in btw))
    o.write(f'# A/B diagnostic positions (mosaic scan): {len(diag)}; NUMT-private positions: {len(priv)}\n')
    o.write('# N placement (observed N at the tested sites vs expected from the sites fraction; one-sided Poisson P)\n')
    for lab, mem, lc, t in ntest:
        nN, nNd, exp, pv = n_enrich(mem, lc, t)
        o.write(f'# {lab}: N {nN}, at tested sites {nNd}, expected {exp:.1f}, enrichment {nNd / exp if exp else float("nan"):.1f}x, P {pv:.1e}\n')

# GenBank annotation check
with open('results/genbank_annotation.tsv', 'w') as o:
    o.write('accession\tCDS_features\tpseudo_flags\tCDS_translations_with_stop\tother_notes\n')
    for rec in open('gb/all.gb').read().split('\n//'):
        m = re.search(r'ACCESSION\s+(\S+)', rec)
        if not m: continue
        ncds = pseudo = stopin = 0; notes = set()
        for f in re.finditer(r'\n     (CDS|gene)\s+(\S+)((?:\n                     .*)*)', rec):
            q = f.group(3)
            if '/pseudo' in q: pseudo += 1
            if f.group(1) != 'CDS': continue
            ncds += 1
            t = re.search(r'/translation="([^"]+)', q, re.S)
            if t and '*' in t.group(1): stopin += 1
            for nn in re.findall(r'/note="([^"]+)', q, re.S):
                nn = ' '.join(nn.split())
                if 'stop codon is completed' not in nn: notes.add(nn[:80])
        o.write(f"{m.group(1)}\t{ncds}\t{pseudo}\t{stopin}\t{'; '.join(sorted(notes))}\n")

# viewer files
order = [n for k in GROUP_ORDER for n in groups.get(k, [])]
lab = {n: f"{n}|{per[n]['group']}|fs{per[n]['fs']}|stop{per[n]['stops']}" for n in order}
# (1) reference-anchored merge of the pairwise alignments, insertion columns kept (frameshifts stay visible).
# ViewAlign marks frameshifts against the column-majority codon phase, so the overview keeps only the 5 B-family copies
# with the best coverage: with all 82 the pseudogene rows outnumber the mitogenomes and the viewer flags their shared
# indels as frameshifts in every real gene. All 82 copies are in viewer_numtB_family_13genes.fasta (with the reference).
best5 = sorted(groups['NUMT_B_family'], key=lambda n: -per[n]['cov'])[:5]
order_full = order
order = [n for n in order_full if per[n]['group'] != 'NUMT_B_family' or n in best5]
def merged_view(order, prefix, per_gene=True):
    rows = {n: [] for n in order}; bed = []; pos = 0
    for g in PCG:
        Lg = len(refseq[g]); w = 0; block = {n: [] for n in order}
        for j in range(Lg + 1):
            mx = max(len(ins[g][n][j]) for n in order if n in ins[g])
            for n in order:
                if n not in proj[g]: block[n].append('-' * (mx + (j < Lg))); continue
                block[n].append(ins[g][n][j].ljust(mx, '-') + (proj[g][n][j] if j < Lg else ''))
            w += mx + (j < Lg)
        with open(f'results/{prefix}_{g}.fasta', 'w') as o:
            for n in order:
                if n in proj[g]: o.write(f'>{lab[n]}\n{"".join(block[n])}\n')
        for n in order: rows[n].append(''.join(block[n]))
        # BED positions are positions in the reference sequence itself: ViewAlign maps them through the backbone row's gaps
        bed.append(f'{lab[REF]}\t{pos}\t{pos + Lg}\t{g}\t0\t+\t{pos}\t{pos + Lg}\t0,0,0\t1\t{Lg}\t0\tCDS {g}')
        pos += Lg
    with open(f'results/{prefix}_13genes.fasta', 'w') as o:
        for n in order: o.write(f'>{lab[n]}\n{"".join(rows[n])}\n')
    open(f'results/{prefix}_13genes.bed', 'w').write('\n'.join(bed) + '\n')

merged_view(order, 'viewer_all')
merged_view([REF] + groups['NUMT_B_family'], 'viewer_numtB_family')
# (2) MACSE multiple alignment of own mtDNA + GenBank mitogenomes ('!' -> '-')
def readfa(p):
    d = collections.OrderedDict(); k = None
    for l in open(p):
        l = l.rstrip('\n')
        if l.startswith('>'): k = l[1:]; d[k] = []
        elif k: d[k].append(l.strip())
    return {a: ''.join(b) for a, b in d.items()}
mord = [n for n in order if n == REF or n.startswith('Sbet_')]
rows = {n: [] for n in mord}; bed = []; pos = 0
for g in PCG:
    m = readfa(f'genes/{g}.main_NT.fna'); w = len(next(iter(m.values())))
    with open(f'results/viewer_mitogenomes_{g}.fasta', 'w') as o:
        for n in mord:
            if n in m: o.write(f'>{lab[n]}\n{m[n].replace("!", "-")}\n')
    for n in mord: rows[n].append(m.get(n, '-' * w).replace('!', '-'))
    L = len(m[REF].replace('-', '').replace('!', ''))
    bed.append(f'{lab[REF]}\t{pos}\t{pos + L}\t{g}\t0\t+\t{pos}\t{pos + L}\t0,0,0\t1\t{L}\t0\tCDS {g}')
    pos += L
with open('results/viewer_mitogenomes_13genes.fasta', 'w') as o:
    for n in mord: o.write(f'>{lab[n]}\n{"".join(rows[n])}\n')
open('results/viewer_mitogenomes_13genes.bed', 'w').write('\n'.join(bed) + '\n')
print(open('results/groups.tsv').read())
print('ANALYZE2_DONE')
