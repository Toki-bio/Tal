#!/usr/bin/env python3
"""Add 'Copies all / firm' and 'Mean identity' columns to the top subfamily table of the six
Talpidae species reports (same columns as SINEderella step6_report.py since commit 115885c).

The pages were built in April by an older step6 and then edited by hand (subfamily letters,
tribe curation), so they are patched in place instead of regenerated.

Data, per species, from the run each page reports on:
  copies:   summary.by_subfam.tsv  total_assigned / firm_assigned
  identity: _sinederella/kit_run_<sp>/plots_data/pctid_<sf>.tsv  (mean of column 2)
Safety: every copy number is checked against the page's own QC table before writing.
Idempotent: a page that already has the columns is left alone.

Run from c:\\work\\Tal:  python _patch_top_table_copies_identity.py
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SPECIES = ['saq', 'ccr', 'toc', 'teu', 'gpy', 'dmo']


def summary_path(sp):
    for p in (ROOT / sp / 'summary.by_subfam.tsv',
              ROOT / '_sinederella' / f'kit_run_{sp}' / 'step2' / 'step2_output' / 'summary.by_subfam.tsv',
              ROOT / '_sinederella' / f'kit_run_{sp}' / 'step2' / 'summary.by_subfam.tsv'):
        if p.is_file():
            return p
    sys.exit(f'{sp}: no summary.by_subfam.tsv')


def load(sp):
    lines = summary_path(sp).read_text(encoding='utf-8').splitlines()
    hdr = lines[0].split('\t')
    i_all, i_firm = hdr.index('total_assigned'), hdr.index('firm_assigned')
    copies = {}
    for l in lines[1:]:
        f = l.split('\t')
        if len(f) > i_all:
            copies[f[0]] = (int(f[i_all]), int(f[i_firm]))
    ident = {}
    for p in (ROOT / '_sinederella' / f'kit_run_{sp}' / 'plots_data').glob('pctid_*.tsv'):
        vals = [float(l.split('\t')[1]) for l in p.read_text(encoding='utf-8').splitlines()
                if l.strip() and len(l.split('\t')) >= 2]
        if vals:
            ident[p.stem[len('pctid_'):]] = (sum(vals) / len(vals), len(vals))
    return copies, ident


def page_qc_copies(html):
    """(firm, total) per subfamily from the page's own copy-number QC table
    (header 'Firm ... Total copies' in newer pages, raw 'firm_assigned ... total_assigned' in older ones)."""
    out = {}
    for m in re.finditer(r"<table.*?</table>", html, flags=re.S):
        rows = re.findall(r'<tr.*?</tr>', m.group(0), flags=re.S)
        if not rows:
            continue
        hdr = [re.sub(r'<[^>]+>', '', c).strip() for c in re.findall(r'<t[hd][^>]*>.*?</t[hd]>', rows[0], flags=re.S)]
        for f_name, t_name in (('Firm', 'Total copies'), ('firm_assigned', 'total_assigned')):
            if f_name in hdr and t_name in hdr:
                i_f, i_t = hdr.index(f_name), hdr.index(t_name)
                for r in rows[1:]:
                    cells = [re.sub(r'<[^>]+>', '', c).strip() for c in re.findall(r'<td[^>]*>.*?</td>', r, flags=re.S)]
                    if len(cells) > max(i_f, i_t):
                        out[cells[0]] = (int(cells[i_f].replace(',', '')), int(cells[i_t].replace(',', '')))
        if out:
            return out
    return out


def patch(sp):
    page = ROOT / sp / 'report.html'
    html = page.read_text(encoding='utf-8')
    if 'Copies all&nbsp;/&nbsp;firm' in html:
        print(f'{sp}: already patched, skipped')
        return
    copies, ident = load(sp)
    qc = page_qc_copies(html)
    if set(qc) != set(copies):
        sys.exit(f'{sp}: page QC table subfamilies {sorted(qc)} != summary file {sorted(copies)}; not patching')
    for sf, (tot, firm) in copies.items():
        if qc[sf] != (firm, tot):
            sys.exit(f'{sp}: {sf} page QC says {qc[sf]} but summary file says {(firm, tot)}; not patching')
    i = html.find('Top 100 by bitscore')
    s = html.rfind('<table', 0, i)
    e = html.find('</table>', i) + len('</table>')
    table = html[s:e]
    th = ("<th title='All = every locus assigned to this subfamily by its best vote (total_assigned); "
          "firm = 10/10 unanimous votes and bitscore above threshold (firm_assigned)'>Copies all&nbsp;/&nbsp;firm</th>"
          "<th title='Mean ssearch36 % identity of copies to the subfamily consensus (copies with 10/10 "
          "unanimous votes, up to 10,000 sampled). Higher = younger.'>Mean identity</th>")
    assert table.count('<th>Subfamily</th>') == 1
    table = table.replace('<th>Subfamily</th>', '<th>Subfamily</th>' + th, 1)

    def add_cells(m):
        row = m.group(0)
        first = re.match(r'<tr><td>(.*?)</td>', row, flags=re.S)
        if not first:
            return row
        names = re.findall(r'\(([A-Za-z0-9_]+)\)', first.group(1)) or re.findall(r'<code>([^<]+)</code>', first.group(1))
        sf = next((n for n in names if n in copies or n in ident), None)
        if sf is None:
            sys.exit(f'{sp}: no data for row {first.group(1)!r}')
        tot, firm = copies[sf]
        mean, n = ident[sf]
        cells = (f"<td style='text-align:right'>{tot:,} / {firm:,}</td>"
                 f"<td style='text-align:right' title='mean of {n:,} copies with unanimous votes'>{mean:.1f}&nbsp;%</td>")
        print(f'  {sp} {sf:<14} {tot:>9,} / {firm:>9,}   {mean:5.1f} %  (n={n:,})')
        return row[:first.end()] + cells + row[first.end():]

    body = re.sub(r'<tr><td>.*?</tr>', add_cells, table, flags=re.S)
    page.write_text(html[:s] + body + html[e:], encoding='utf-8', newline='')
    print(f'{sp}: patched')


if __name__ == '__main__':
    for sp in (sys.argv[1:] or SPECIES):
        patch(sp)
