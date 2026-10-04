#!/usr/bin/env python3
"""GenBank flat file -> BED12+1 for the ViewAlign annotation track.

    python3 gb2bed.py record.gb [chrom_name] > genes.bed

chrom_name is the sequence name as it appears in the alignment (default: the
LOCUS name). Features kept: CDS, tRNA, rRNA, rep_origin, D-loop and
misc_feature (control region). `gene` features are skipped (duplicates of the
CDS/tRNA/rRNA they wrap). Columns 1-12 are standard BED (0-based, half-open,
itemRgb by feature class); column 13 is a free-text description shown in the
track tooltip (feature type and /product).
"""
import re
import sys

COLORS = {               # itemRgb per feature class
    'tRNA': '124,179,66',
    'rRNA': '251,140,0',
    'ND': '30,136,229',
    'COX': '0,151,167',
    'ATP': '126,87,194',
    'CYTB': '57,73,171',
    'CR': '158,158,158',
    'OL': '189,189,189',
    'other': '120,144,156',
}


def parse_features(path):
    feats, cur = [], None
    in_feat = False
    for line in open(path, encoding='utf-8', errors='replace'):
        if line.startswith('FEATURES'):
            in_feat = True
            continue
        if line.startswith('ORIGIN') or line.startswith('//'):
            break
        if not in_feat:
            continue
        m = re.match(r'^     (\S+)\s+(\S.*)$', line)
        if m:
            cur = {'type': m.group(1), 'loc': m.group(2).strip(), 'q': {}}
            feats.append(cur)
            continue
        m = re.match(r'^\s{21}/(\w+)(?:=("?)(.*?)\2)?\s*$', line)
        if m and cur is not None:
            cur['q'][m.group(1)] = m.group(3) if m.group(3) is not None else ''
            cur['_last'] = m.group(1)
        elif cur is not None and cur.get('_last') and line.startswith(' ' * 21):
            # continuation of a multi-line qualifier (e.g. a long /note)
            cur['q'][cur['_last']] += ' ' + line.strip().rstrip('"')
        elif cur is not None and line.startswith(' ' * 21):
            cur['loc'] += line.strip()     # continuation of a long location
    return feats


def span(loc):
    strand = '-' if loc.startswith('complement(') else '+'
    nums = [int(n) for n in re.findall(r'\d+', loc)]
    return min(nums), max(nums), strand


def name_and_class(f):
    t, q = f['type'], f['q']
    gene, product = q.get('gene', ''), q.get('product', '')
    if t == 'CDS':
        g = (gene or product).upper().replace('COB', 'CYTB')
        for k in ('ND', 'COX', 'ATP', 'CYTB'):
            if g.startswith(k):
                return gene or product, k
        return gene or product, 'other'
    if t == 'tRNA':
        return product or gene, 'tRNA'
    if t == 'rRNA':
        p = product or gene
        short = re.sub(r'\s*ribosomal RNA', ' rRNA', p)
        return short, 'rRNA'
    if t == 'rep_origin':
        return 'OL', 'OL'
    if t in ('D-loop', 'misc_feature'):
        note = (q.get('note', '') + ' ' + product).lower()
        if 'control' in note or 'd-loop' in note or t == 'D-loop':
            return 'control region', 'CR'
        return None, None
    return None, None


def main():
    path = sys.argv[1]
    locus = re.match(r'LOCUS\s+(\S+)', open(path, encoding='utf-8', errors='replace').readline()).group(1)
    chrom = sys.argv[2] if len(sys.argv) > 2 else locus
    print('track name="%s genes" description="GenBank annotation of %s" itemRgb="On"' % (locus, locus))
    for f in parse_features(path):
        name, cls = name_and_class(f)
        if not name:
            continue
        a, b, strand = span(f['loc'])
        desc = f['type'] + (': ' + f['q']['product'] if f['q'].get('product') else '')
        if f['q'].get('note') and f['type'] not in ('CDS',):
            desc += '; ' + f['q']['note']
        print('\t'.join(str(x) for x in [chrom, a - 1, b, name, 0, strand, a - 1, b, COLORS[cls],
                                          1, b - a + 1, 0, desc]))


if __name__ == '__main__':
    main()
