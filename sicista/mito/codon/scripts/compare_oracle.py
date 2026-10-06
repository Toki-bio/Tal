#!/usr/bin/env python3
# ViewAlign MACSE port vs real MACSE v2.07 on the same gene inputs (gc 2, one reliable + less-reliable sequences).
# Identity = same names, same row order, every '!' and '-' the same; also the AA output.
import os, collections
os.chdir(os.path.expanduser('~/Sicista2026/mito/pseudo_check'))


def readfa(p):
    d = collections.OrderedDict(); n = None
    for l in open(p):
        l = l.rstrip('\n')
        if l.startswith('>'): n = l[1:]; d[n] = []
        elif n is not None: d[n].append(l.strip())
    return collections.OrderedDict((k, ''.join(v)) for k, v in d.items())


same = 0; n = 0
for g in ['ATP8', 'ND4L', 'ND3', 'ND6', 'COX2', 'ATP6', 'COX3', 'ND1', 'ND2', 'CYTB', 'ND4', 'COX1', 'ND5']:
    pm = f'oracle/{g}.macse_NT.fna'
    if not os.path.exists(pm) or not os.path.exists(f'genes/{g}.port_NT.fna'):
        print(g, 'not yet'); continue
    a = readfa(f'genes/{g}.port_NT.fna'); b = readfa(pm)
    aa_a = readfa(f'genes/{g}.port_AA.faa'); aa_b = readfa(f'oracle/{g}.macse_AA.faa')
    order = list(a) == list(b)
    diff = [k for k in b if a.get(k) != b[k]]
    diffaa = [k for k in aa_b if aa_a.get(k) != aa_b[k]]
    ok = order and not diff and not diffaa and set(a) == set(b)
    n += 1; same += ok
    print(f'{g}\t{"IDENTICAL" if ok else "DIFFERENT"}\trows {len(a)}/{len(b)}\twidth {len(next(iter(a.values())))}/{len(next(iter(b.values())))}'
          f'\torder {"same" if order else "differs"}\tNT diffs {len(diff)}\tAA diffs {len(diffaa)}\t{",".join(diff[:3])}')
print(f'identical {same}/{n}')
