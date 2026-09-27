#!/usr/bin/env python3
"""Build chiroptera.html, the Chiroptera SINE project page, grouped by clade.

Inputs (all in the Tal repo):
  chiroptera/genomes.tsv        code, species, family, accession, assembly, suborder, superfamily,
                                prior_evidence; rows in the order of Hao et al. 2023 Fig. 3, top to
                                bottom. code "-" = a family with no genome searched (a note, no card).
  <code>/summary.by_subfam.tsv  SINEderella step3 summary. Rhin-1 and VES are the original bank.
                                MEG-RL, MEG-RS, MEG-T2, MEG-TR appear after the Gogolevsky 2009 add;
                                a summary with none of those four names means MEG was not searched.
  <code>/alignments/*.aln.fa    published alignments (SubFam input, de novo candidates)
  LEGACY below                  the four genomes searched before the one-per-family round, whose
                                published files do not follow that layout

The whole page is regenerated; run this after every publish instead of editing the HTML by hand.
"""
import csv
import glob
import html
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAGE = os.path.join(ROOT, "chiroptera.html")
SINES = ["Rhin-1", "VES"]
# Gogolevsky, Vassetzky & Kramerov 2009, Genomics 93:494-500. Four families, not one consensus.
# The tree badge is their sum: assignment gives each locus to one family, so the sum does not double-count.
MEG = ["MEG-RL", "MEG-RS", "MEG-T2", "MEG-TR"]
MEG_SHORT = {"MEG-RL": "RL", "MEG-RS": "RS", "MEG-T2": "T2", "MEG-TR": "TR"}
DENOVO_PENDING = set()  # hla, tbr de novo chains finished 2026-09-27

# (copies, firm, sim median) per SINE; None = not searched. Earlier runs: Rhin-1 was the
# 177 bp sanitised consensus (IUPAC codes deleted), and rsi/rre/rda had no VES in the bank.
LEGACY = {
    "rsi": {"Rhin-1": (30079, 23545, 0.48), "VES": None, "tag": ("tag-ok", "full report"),
            "note": "Rhin-1-only search; peel bank r1&ndash;r10: 10 subfamilies &middot; 67,162 assigned copies",
            "acts": [("btn", "rsi/report.html", "<i>Rhinolophus sinicus</i> &mdash; full report"),
                     ("msa", "rsi/alignments/rsi_consensuses.aln.fa", "rsi all consensi (r1-r10 + Rhin-1)", "All Consensi in MSA"),
                     ("msa", "rhin/alignments/rsi_subfam_input_30k.aln.fa", "rsi Rhin-1 SubFam input, 600 chunk consensi + anchor", "SubFam input (600 chunks)"),
                     ("sec", "rhin.html", "Peel alignments")]},
    "rre": {"Rhin-1": (32410, 23073, 0.48), "VES": None, "tag": ("tag-ok", "full report"),
            "note": "Rhin-1 + MEG search (no VES)",
            # full report with alignments since the MEG add (2026-09-27)
            "acts": [("btn", "rre/report.html", "<i>Rhinolophus rex</i> &mdash; full report"),
                     ("msa", "rhin/alignments/rre_subfam_input_30k.aln.fa", "rre Rhin-1 SubFam input, 600 chunk consensi + anchor", "SubFam input (600 chunks)")]},
    "rda": {"Rhin-1": (30169, 23850, 0.48), "VES": None, "tag": ("tag-partial", "pending manual review"),
            "note": "Rhin-1-only search",
            "acts": [("msa", "rhin/alignments/rda_subfam_input_30k.aln.fa", "rda (hap1) Rhin-1 SubFam input, 600 chunk consensi + anchor", "SubFam input (600 chunks)")]},
    "lly": {"Rhin-1": (290, 290, 0.16), "VES": (11, 11, 0.14), "tag": ("tag-ok", "full report"),
            "note": "fragmented assembly (1.9 M scaffolds); <i>Macroderma</i> (mgi) is the better Megadermatidae genome",
            # full report with alignments since the MEG add (2026-09-27); the old chiroptera/lly report had none
            "acts": [("btn", "lly/report.html", "<i>Lyroderma lyra</i> &mdash; full report"),
                     ("msa", "chiroptera/lly/lly_subfam_input.aln.fa", "lly SubFam input, 6 chunk consensi + Rhin-1 + VES", "SubFam input (6 chunks)")]},
}


def fmt(n):
    return f"{int(n):,}"


def esc(s):
    return html.escape(s, quote=True)


def msa_btn(path, title, label, cls="secondary"):
    return (f'<a class="btn {cls}" href="javascript:void(0)" '
            f"onclick=\"openMSA('{path}','{esc(title)}')\">{label}</a>")


def read_summary(code):
    p = os.path.join(ROOT, code, "summary.by_subfam.tsv")
    if not os.path.exists(p):
        return None
    with open(p, encoding="utf8") as fh:
        return {r["subfam"]: r for r in csv.DictReader(fh, delimiter="\t")}


def sine_tuple(rs, name):
    r = rs[name]
    return (int(r["total_assigned"]), int(r["firm_assigned"]), float(r["sim_median"]))


def meg_searched(vals):
    return any(s in vals for s in MEG)


def meg_sum(vals):
    if not meg_searched(vals):
        return None
    return sum(vals[s][0] for s in MEG if vals.get(s))


def load(g):
    """Return (status, {sine: (copies, firm, sim) or None}, tag, note, actions).

    MEG keys are present only once that genome's summary contains a MEG family.
    Absent keys mean not searched, which is different from a searched family with 0 copies.
    """
    c, sp = g["code"], g["species"]
    if c in LEGACY:
        L = LEGACY[c]
        acts = []
        for a in L["acts"]:
            if a[0] == "btn":
                acts.append(f'<a class="btn" href="{a[1]}">{a[2]}</a>')
            elif a[0] == "sec":
                acts.append(f'<a class="btn secondary" href="{a[1]}">{a[2]}</a>')
            else:
                acts.append(msa_btn(a[1], a[2], a[3]))
        vals = {s: L[s] for s in SINES}
        rs = read_summary(c)
        if rs and any(s in rs for s in MEG):
            for s in SINES:
                if s in rs:
                    vals[s] = sine_tuple(rs, s)
            for s in MEG:
                vals[s] = sine_tuple(rs, s) if s in rs else (0, 0, None)
        return "done", vals, L["tag"], L["note"], acts
    rs = read_summary(c)
    if rs is None:
        return "running", {}, ("tag-no", "running"), "", []
    vals = {s: (sine_tuple(rs, s) if s in rs else (0, 0, None)) for s in SINES}
    if any(s in rs for s in MEG):
        for s in MEG:
            vals[s] = sine_tuple(rs, s) if s in rs else (0, 0, None)
    acts = [f'<a class="btn" href="{c}/report.html"><i>{sp}</i> &mdash; full report</a>']
    sub = f"{c}/alignments/{c}_subfam_input_30k.aln.fa"
    if os.path.exists(os.path.join(ROOT, sub)):
        n = sum(1 for l in open(os.path.join(ROOT, sub), encoding="utf8") if l.startswith(">")) - 2
        acts.append(msa_btn(sub, f"{c} ({sp}) SubFam input, {n} chunk consensi + Rhin-1 + VES",
                            f"SubFam input ({n} chunk{'s' if n != 1 else ''})"))
    for dn in sorted(glob.glob(os.path.join(ROOT, c, "alignments", f"{c}_denovo_candidates_*chunks.aln.fa"))):
        k = re.search(r"_(\d+)chunks", dn).group(1)
        acts.append(msa_btn(f"{c}/alignments/{os.path.basename(dn)}",
                            f"{c} de novo candidates (AnnoSINE2 + fragment scan), {k} chunk consensi - for manual review",
                            f"De novo candidates ({k} chunks)", "green"))
    return "done", vals, ("tag-ok", "full report"), "", acts


def tal_style():
    """The Tal main page's own <style> block, copied verbatim so both pages share one design."""
    idx = open(os.path.join(ROOT, "index.html"), encoding="utf8").read()
    return idx[idx.index("<style>"):idx.index("</style>") + len("</style>")]


def copies_tag(v, status):
    if status == "running":
        return '<span class="tag tag-no">running</span>'
    if v is None:
        return '<span class="tag tag-no">not searched</span>'
    if v[0] == 0:
        return '<span class="tag tag-no">0</span>'
    cls = "tag-ok" if v[0] >= 1000 else "tag-partial"
    return f'<span class="tag {cls}">{fmt(v[0])} &middot; sim {v[2]:.2f}</span>'


def meg_cell(vals, status):
    """One table cell for the four MEG families. Sum colours the tag; the text lists each family."""
    if status == "running":
        return '<span class="tag tag-no">running</span>'
    if not meg_searched(vals):
        return '<span class="tag tag-no">not searched</span>'
    bits = []
    for s in MEG:
        v = vals.get(s)
        label = MEG_SHORT[s]
        if v is None or v[0] == 0:
            bits.append(f"{label} 0")
        else:
            bits.append(f"{label} {fmt(v[0])} sim {v[2]:.2f}")
    n = meg_sum(vals)
    cls = "tag-ok" if n >= 1000 else "tag-partial" if n else "tag-no"
    return f'<span class="tag {cls}">{" &middot; ".join(bits)}</span>'


# Hao et al. 2023 Fig. 3 topology. Internal node = (age Ma, [children]); leaf = (family, crown age or None).
TREE = (61.4, [
    (54.8, [("Pteropodidae", 28.5),
            (48.4, [(38.3, [(35.2, [("Rhinonycteridae", None), ("Rhinolophidae", 17.7)]),
                            ("Hipposideridae", 14.7)]),
                    (45.8, [(41.0, [("Megadermatidae", 20.2), ("Craseonycteridae", None)]),
                            ("Rhinopomatidae", None)])])]),
    (58.9, [(57.2, [("Myzopodidae", None),
                    (54.7, [("Nycteridae", 15.4), ("Emballonuridae", 46.7)])]),
            (58.0, [(50.6, [(47.2, [(43.0, [("Phyllostomidae", 36.8), ("Mormoopidae", 21.0)]),
                                    (44.0, [(31.1, [("Noctilionidae", 3.7), ("Furipteridae", None)]),
                                            ("Thyropteridae", None)])]),
                            ("Mystacinidae", None)]),
                    (54.9, [(53.0, [(51.1, [(46.5, [("Cistugidae", None), ("Vespertilionidae", 39.0)]),
                                            ("Miniopteridae", 10.3)]),
                                    ("Molossidae", 24.8)]),
                            ("Natalidae", 11.8)])])])])
RED_DOTS = {54.8, 58.9}  # MRCAs of the two suborders, as in the figure
BANDS = {"Pteropodoidea": "#f3ea7a", "Rhinolophoidea": "#a8d69c", "Emballonuroidea": "#f7c29a",
         "Noctilionoidea": "#a9dfe3", "Vespertilionoidea": "#d9b3d9"}
SUBORDER = {"Yinpterochiroptera": "#2b7a3b", "Yangochiroptera": "#1c3c6e"}
SINE_COL = {"Rhin-1": "#0070C0", "VES": "#E00000", "MEG": "#6A3D9A"}  # Rhin and VES match the annotated figure; MEG does not appear on it
SINE_PALE = {"Rhin-1": "#d7e8f7", "VES": "#fbd6d6", "MEG": "#EDE4F7"}  # opaque tints, so clade bands do not show through
EPOCHS = [  # (row, name, from Ma, to Ma, ICS colour)
    ("Epoch", "Paleocene", 66, 56, "#FDA75F"), ("Epoch", "Eocene", 56, 33.9, "#FDB46C"),
    ("Epoch", "Oligocene", 33.9, 23.03, "#FDC07A"), ("Epoch", "Miocene", 23.03, 5.33, "#FFFF00"),
    ("Epoch", "Pl.", 5.33, 2.58, "#FFFF99"), ("Epoch", "", 2.58, 0, "#FFF2AE"),
    ("Period", "Paleogene", 66, 23.03, "#FD9A52"), ("Period", "Neogene", 23.03, 2.58, "#FFE619"),
    ("Period", "Qu.", 2.58, 0, "#F9F97F"), ("Era", "Cenozoic", 66, 0, "#F2F91D")]


def short(n):
    return f"{n / 1000:.0f}k" if n >= 10000 else f"{n / 1000:.1f}k" if n >= 1000 else str(n)


def report_href(g, status):
    if g["code"] in LEGACY:
        return next((a[1] for a in LEGACY[g["code"]]["acts"] if a[0] == "btn"), None)
    return f'{g["code"]}/report.html' if status == "done" else None


def tree_svg(rows):
    """Hao et al. Fig. 3 redrawn; each leaf carries this project's Rhin-1 / VES result."""
    leaves = []
    for g in rows:
        if g["family"] not in leaves:
            leaves.append(g["family"])
    STEP, TOP, X0, S = 34, 22, 18, 8.0
    x = lambda age: round(X0 + (66 - age) * S, 1)
    tipx = x(0)
    ly = {f: TOP + STEP / 2 + i * STEP for i, f in enumerate(leaves)}
    H_TREE = TOP + STEP * len(leaves)
    out, defs = [], []

    # clade bands, fading in from the left as in the figure, and the suborder bars
    for sf, col in BANDS.items():
        fams = [g["family"] for g in rows if g["superfamily"] == sf]
        y0, y1 = ly[fams[0]] - STEP / 2, ly[fams[-1]] + STEP / 2
        defs.append(f'<linearGradient id="band_{sf}" x1="0" x2="1"><stop offset="0" stop-color="{col}" stop-opacity="0"/>'
                    f'<stop offset=".45" stop-color="{col}" stop-opacity=".35"/><stop offset="1" stop-color="{col}"/></linearGradient>')
        out.append(f'<rect x="{X0}" y="{y0}" width="{1212 - X0}" height="{y1 - y0}" fill="url(#band_{sf})"/>')
        if y1 - y0 > 60:
            fs = min(15, round((y1 - y0 - 8) / (0.56 * len(sf)), 1))  # shrink to fit short bands
            out.append(f'<text transform="translate(1198,{(y0 + y1) / 2}) rotate(90)" text-anchor="middle" '
                       f'dominant-baseline="central" font-size="{fs}" font-weight="700" fill="#222">{sf}</text>')
    for so, col in SUBORDER.items():
        fams = [g["family"] for g in rows if g["suborder"] == so]
        y0, y1 = ly[fams[0]] - STEP / 2, ly[fams[-1]] + STEP / 2
        out.append(f'<rect x="1214" y="{y0 + 1}" width="42" height="{y1 - y0 - 2}" rx="12" fill="{col}"/>')
        out.append(f'<text transform="translate(1235,{(y0 + y1) / 2}) rotate(90)" text-anchor="middle" '
                   f'dominant-baseline="central" font-size="19" font-weight="700" fill="#fff">{so}</text>')

    # branches, crown triangles, node ages, red dots
    def draw(node):
        if isinstance(node[0], str):
            fam, crown = node
            y = ly[fam]
            if crown:
                out.append(f'<polygon points="{x(crown)},{y} {tipx},{y - 12} {tipx},{y + 12}" '
                           f'fill="#fff" stroke="#000" stroke-width="1.2"/>')
                return x(crown), y
            return tipx, y
        age, kids = node
        pts = [draw(k) for k in kids]
        nx, ys = x(age), [p[1] for p in pts]
        for kx, ky in pts:
            out.append(f'<line x1="{nx}" y1="{ky}" x2="{kx}" y2="{ky}" stroke="#000" stroke-width="1.4"/>')
        out.append(f'<line x1="{nx}" y1="{min(ys)}" x2="{nx}" y2="{max(ys)}" stroke="#000" stroke-width="1.4"/>')
        ny = (min(ys) + max(ys)) / 2
        out.append(f'<text x="{nx + 3}" y="{ny - 4}" font-size="10.5" font-weight="700" fill="#222">{age}</text>')
        if age in RED_DOTS:
            out.append(f'<circle cx="{nx}" cy="{ny}" r="4.5" fill="#e8202a"/>')
        return nx, ny
    rx, ry = draw(TREE)
    out.append(f'<line x1="{x(65)}" y1="{ry}" x2="{rx}" y2="{ry}" stroke="#000" stroke-width="1.4"/>')

    # leaves: family name, one badge per SINE, genome codes linking to reports.
    # MEG is one badge: the sum of the four families. Each locus is assigned to one family.
    BX = {"Rhin-1": 992, "VES": 1056, "MEG": 1120}
    LABEL = {"Rhin-1": "Rhin", "VES": "VES", "MEG": "MEG"}
    for s, bx in BX.items():
        out.append(f'<text x="{bx + 30}" y="{TOP - 6}" text-anchor="middle" font-size="12" font-weight="700" '
                   f'fill="{SINE_COL[s]}">{LABEL[s]}</text>')
    for fam in leaves:
        y = ly[fam]
        gs = [g for g in rows if g["family"] == fam and g["code"] != "-"]
        loaded = [(g, load(g)) for g in gs]
        out.append(f'<text x="{tipx + 10}" y="{y}" dominant-baseline="central" font-size="16" fill="#111">{fam}</text>')
        for s, bx in BX.items():
            col = SINE_COL[s]
            if s == "MEG":
                vals = []
                for g, r in loaded:
                    if r[0] == "done" and meg_searched(r[1]):
                        parts = []
                        for m in MEG:
                            v = r[1].get(m)
                            if v and v[0]:
                                parts.append(f"{MEG_SHORT[m]}: {fmt(v[0])} copies, sim {v[2]:.2f}")
                            else:
                                parts.append(f"{MEG_SHORT[m]}: 0")
                        vals.append((g["code"], meg_sum(r[1]), "; ".join(parts)))
                running = any(r[0] == "running" for g, r in loaded)
                tip = "; ".join(f"{c}: {detail}" for c, _n, detail in vals)
                if vals:
                    n = max(n for _, n, _ in vals)
                    if n >= 1000:
                        fill, stroke, tc, label = col, col, "#fff", short(n)
                    elif n > 0:
                        fill, stroke, tc, label = SINE_PALE[s], col, col, short(n)
                    else:
                        fill, stroke, tc, label = "#fff", "#999", "#777", "0"
                elif running:
                    fill, stroke, tc, label, tip = "#f4f4f4", "#ccc", "#999", "…", "search running"
                else:
                    fill, stroke, tc, label, tip = "#e4e4e4", "#bbb", "#777", "NA", "not searched"
            else:
                vals = [(g["code"], r[1][s]) for g, r in loaded if r[0] == "done" and r[1].get(s) is not None]
                running = any(r[0] == "running" for g, r in loaded)
                tip = "; ".join(f"{c}: {fmt(v[0])} copies" + (f", sim {v[2]:.2f}" if v[0] else "") for c, v in vals)
                if vals:
                    n = max(v[0] for _, v in vals)
                    if n >= 1000:
                        fill, stroke, tc, label = col, col, "#fff", short(n)
                    elif n > 0:
                        fill, stroke, tc, label = SINE_PALE[s], col, col, short(n)
                    else:
                        fill, stroke, tc, label = "#fff", "#999", "#777", "0"
                elif running:
                    fill, stroke, tc, label, tip = "#f4f4f4", "#ccc", "#999", "…", "search running"
                else:
                    fill, stroke, tc, label, tip = "#e4e4e4", "#bbb", "#777", "NA", "not searched"
            out.append(f'<g><title>{fam} {s}: {tip}</title>'
                       f'<rect x="{bx}" y="{y - 10}" width="60" height="20" rx="4" fill="{fill}" stroke="{stroke}"/>'
                       f'<text x="{bx + 30}" y="{y}" text-anchor="middle" dominant-baseline="central" '
                       f'font-size="12" font-weight="700" fill="{tc}">{label}</text></g>')
        # species names, genus abbreviated after the first of the same genus; each links to its report
        cx, prev_genus = 720, None
        for g, r in loaded:
            href = report_href(g, r[0])
            genus, _, epi = g["species"].partition(" ")
            name = f"{genus[0]}. {epi}" if genus == prev_genus else g["species"]
            prev_genus = genus
            t = (f'<text x="{cx}" y="{y}" dominant-baseline="central" font-size="13" font-style="italic" '
                 f'fill="{"#4C72B0" if href else "#555"}"{" text-decoration=" + chr(34) + "underline" + chr(34) if href else ""}>'
                 f'<title>{g["code"]}</title>{name}</text>')
            out.append(f'<a href="{href}">{t}</a>' if href else t)
            cx += round(6.2 * len(name)) + 14

    # geological time bars and axis
    ya = H_TREE + 8
    rowy = {"Epoch": ya, "Period": ya + 16, "Era": ya + 32}
    for row, name, a, b, col in EPOCHS:
        out.append(f'<rect x="{x(a)}" y="{rowy[row]}" width="{round(x(b) - x(a), 1)}" height="16" '
                   f'fill="{col}" stroke="#fff" stroke-width=".6"/>')
        if name:
            out.append(f'<text x="{(x(a) + x(b)) / 2}" y="{rowy[row] + 8}" text-anchor="middle" '
                       f'dominant-baseline="central" font-size="10.5">{name}</text>')
    for row, yy in rowy.items():
        out.append(f'<text x="{tipx + 8}" y="{yy + 8}" dominant-baseline="central" font-size="10.5">{row}</text>')
    yt = ya + 50
    out.append(f'<line x1="{x(66)}" y1="{yt}" x2="{tipx}" y2="{yt}" stroke="#000"/>')
    for t in [66, 60, 50, 40, 30, 20, 10, 0]:
        out.append(f'<line x1="{x(t)}" y1="{yt}" x2="{x(t)}" y2="{yt + 5}" stroke="#000"/>'
                   f'<text x="{x(t)}" y="{yt + 17}" text-anchor="middle" font-size="11">{t}</text>')
    out.append(f'<text x="{tipx + 8}" y="{yt + 17}" font-size="11">(Ma)</text>')
    H = yt + 26
    return (f'<svg viewBox="0 0 1280 {H}" width="100%" style="min-width:1040px;display:block;" '
            f'xmlns="http://www.w3.org/2000/svg" role="img" '
            f'aria-label="Bat family tree (Hao et al. 2023) with Rhin-1, VES and MEG copy numbers per family">'
            f'<defs>{"".join(defs)}</defs>{"".join(out)}</svg>')


NAV_START, NAV_END = "<!-- BACKNAV START -->", "<!-- BACKNAV END -->"


def add_back_nav(report_rel):
    """Put Chiroptera / Tal main page buttons at the top of a species report.

    Reports are regenerated by SINEderella on the server, so this runs on every build;
    the marked block is replaced, never duplicated."""
    path = os.path.join(ROOT, report_rel)
    if not os.path.exists(path):
        return False
    up = "../" * report_rel.count("/")
    btn = ("display:inline-block;background:{bg};color:#fff;padding:4px 12px;border-radius:4px;"
           "font-size:.84rem;text-decoration:none;margin:0 6px 10px 0;")
    nav = (f'{NAV_START}\n  <div>'
           f'<a href="{up}chiroptera.html" style="{btn.format(bg="#2d2d8f")}">&larr; Chiroptera SINEs</a>'
           f'<a href="{up}index.html" style="{btn.format(bg="#4C72B0")}">Tal SINE main page</a>'
           f'</div>\n  {NAV_END}')
    page = open(path, encoding="utf8").read()
    if NAV_START in page:
        new = re.sub(re.escape(NAV_START) + ".*?" + re.escape(NAV_END), lambda m: nav, page, flags=re.S)
    elif "<header>" in page:
        new = page.replace("<header>", "<header>\n  " + nav, 1)
    else:
        new = page.replace("<body>", "<body>\n" + nav, 1)
    if new != page:
        open(path, "w", encoding="utf8", newline="\n").write(new)
    return True


def main():
    with open(os.path.join(ROOT, "chiroptera", "genomes.tsv"), encoding="utf8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    for g in rows:
        if g["code"] == "-":
            continue
        rep = (next((a[1] for a in LEGACY[g["code"]]["acts"] if a[0] == "btn"), None)
               if g["code"] in LEGACY else f'{g["code"]}/report.html')
        if rep:
            add_back_nav(rep)

    groups = []  # (suborder, superfamily) in tree order
    for g in rows:
        k = (g["suborder"], g["superfamily"])
        if k not in groups:
            groups.append(k)

    sections, trs = [], []
    n_done = n_gen = 0
    for so, sf in groups:
        cards, notes = [], []
        first = True
        for g in rows:
            if g["superfamily"] != sf:
                continue
            c, sp, fam = g["code"], g["species"], g["family"]
            if c == "-":
                notes.append(f"{fam} is not searched here: VES was described from it, so it was left out.")
                continue
            border = ' style="border-top:3px solid #2d2d8f;"' if first else ""
            first = False
            n_gen += 1
            status, vals, tag, note, acts = load(g)
            n_done += status == "done"
            prior = esc(g["prior_evidence"])

            meta = []
            if status == "done":
                bits = [f"{s} {fmt(v[0])}" for s, v in vals.items() if s in SINES and v is not None]
                if meg_searched(vals):
                    bits.append(f"MEG {fmt(meg_sum(vals))}")
                meta.append(" &middot; ".join(bits) + " copies")
            meta.append(f'<a href="https://www.ncbi.nlm.nih.gov/datasets/genome/{g["accession"]}/" '
                        f'target="_blank">{g["accession"]}</a>')
            extra = ([f"prior evidence: {prior}"] if prior else []) + ([note] if note else [])
            meg_line = ""
            if meg_searched(vals):
                meg_line = " &middot; ".join(
                    f"{MEG_SHORT[s]} {fmt(vals[s][0])}" + (f" sim {vals[s][2]:.2f}" if vals[s][0] else "")
                    for s in MEG)
            card = (f'    <div class="sp-card">\n'
                    f'      <h3>{c} <span class="tag {tag[0]}">{tag[1]}</span></h3>\n'
                    f'      <div class="sci">{sp} &mdash; {fam}</div>\n'
                    f'      <div style="font-size:.82rem;color:var(--muted);margin-top:4px;">\n'
                    f'        {" &middot; ".join(meta)}\n'
                    f'      </div>\n')
            if extra:
                card += (f'      <div style="font-size:.82rem;color:var(--muted);margin-top:4px;">'
                         f'{"; ".join(extra)}</div>\n')
            if meg_line:
                card += (f'      <div style="font-size:.82rem;color:var(--muted);margin-top:4px;">'
                         f'MEG {meg_line}</div>\n')
            if acts:
                card += '      <div class="actions">\n        ' + "\n        ".join(acts) + '\n      </div>\n'
            cards.append(card + '    </div>\n')

            href = (next((a[1] for a in LEGACY[c]["acts"] if a[0] == "btn"), None) if c in LEGACY
                    else f"{c}/report.html")
            rep = (f'<span class="tag tag-ok"><a href="{href}">HTML</a></span>' if status == "done" and href
                   else '<span class="tag tag-no">&mdash;</span>')
            if status != "done":
                aln = '<span class="tag tag-no">&mdash;</span>'
            elif c in ("rre", "rda", "lly"):
                aln = '<span class="tag tag-partial">SubFam input</span>'
            else:
                aln = '<span class="tag tag-ok">top100 · rand100 · subfam</span>'
            if glob.glob(os.path.join(ROOT, c, "alignments", f"{c}_denovo_candidates_*chunks.aln.fa")):
                dn = '<span class="tag tag-ok">ready</span>'
            elif c in DENOVO_PENDING:
                dn = '<span class="tag tag-partial">running</span>'
            else:
                dn = '<span class="tag tag-no">&mdash;</span>'
            search = ('<span class="tag tag-ok">✓</span>' if status == "done"
                      else '<span class="tag tag-no">running</span>')
            trs.append(
                f'      <tr{border}>\n'
                f'        <td><strong>{c}</strong> <em>{sp}</em> '
                f'<span class="tag tag-other-family" style="margin-left:4px;">{fam}</span></td>\n'
                f'        <td>{prior or "&mdash;"}</td>\n'
                f'        <td>{search}</td>\n'
                f'        <td>{copies_tag(vals.get("Rhin-1"), status)}</td>\n'
                f'        <td>{copies_tag(vals.get("VES"), status)}</td>\n'
                f'        <td>{meg_cell(vals, status)}</td>\n'
                f'        <td>{rep}</td>\n'
                f'        <td>{aln}</td>\n'
                f'        <td>{dn}</td>\n'
                f'      </tr>')
        note_html = "".join(f'  <p style="color:var(--muted);font-size:.9rem;margin-top:4px;">{n}</p>\n'
                            for n in notes)
        sections.append(
            f'<section class="card">\n'
            f'  <h2>{sf} &mdash; {so}</h2>\n'
            f'{note_html}'
            f'  <div class="species-grid">\n\n' + "\n".join(cards) + '\n  </div>\n</section>\n')

    page = f"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Chiroptera SINE — bat analysis</title>
{tal_style()}
</head><body>
<header>
  <h1>Chiroptera SINE — bat analysis</h1>
  <div class="sub">Rhin-1, VES and the four MEG families (SINEbase) searched in one genome per bat family, plus three <i>Rhinolophus</i> genomes.
  Species are grouped by clade following the family tree of Hao et al. 2023 (<i>Integrative Zoology</i> 19:989&ndash;998, Fig. 3).
  <a href="index.html" style="color:#fff;">&larr; Tal SINE main page</a></div>
</header>
<main>

<section class="card">
  <h2>Family tree &mdash; Rhin-1, VES and MEG per family</h2>
  <p style="margin:0 0 10px;font-size:.9rem;color:var(--muted);">Topology, node ages and colours redrawn from Hao et al. 2023, Fig. 3
  (MrBayes time tree; triangles = collapsed families from their crown age). Each leaf shows the copies found in that family's genome(s):
  <b style="color:#0070C0;">Rhin</b> = Rhin-1, <b style="color:#E00000;">VES</b> = VES,
  <b style="color:#6A3D9A;">MEG</b> = MEG-RL + MEG-RS + MEG-T2 + MEG-TR (Gogolevsky et al. 2009);
  solid = 1,000 copies or more, pale = 1&ndash;999,
  NA = not searched, &hellip; = search running. Where a family has several genomes the largest count is shown (hover for all).
  The MEG badge is the sum of the four families: each locus is assigned to one family, so a copy is not counted twice.
  Similarity is printed with each count. A high count at low similarity is a weak match, not a presence call.
  MEG-RS is nearly the 5&prime; part of MEG-RL, so the report's conflict column, not the badge, is where those two overlap.
  Species names link to their reports. Bank: SINEbase Rhin-1 (182 bp) and VES (220 bp), IUPAC codes resolved to a base;
  the earlier <i>Rhinolophus</i> runs used Rhin-1 alone (177 bp, codes deleted).
  MEG: SINEbase MEG-RL (207 bp), MEG-RS (135 bp), MEG-T2 (232 bp), MEG-TR (196 bp);
  MEG-T2 K/R resolved to T/A/A as in the paper's tandem linker and AT-rich tail.
  MEG was added with <code>--add</code>, without searching Rhin-1 and VES again.</p>
  <div style="overflow-x:auto;">
{tree_svg(rows)}
  </div>
</section>

{chr(10).join(sections)}
<section class="card">
  <h2>Analysis Status</h2>
  <table class="tbl">
    <thead>
      <tr><th>Species</th><th>Prior evidence</th><th>Search</th><th>Rhin-1 &middot; copies</th><th>VES &middot; copies</th><th>MEG &middot; copies</th><th>Report</th><th>Alignments</th><th>De novo</th></tr>
    </thead>
    <tbody>
{chr(10).join(trs)}
    </tbody>
  </table>
  <p style="color:var(--muted);font-size:.82rem;margin-top:8px;">
    {n_done} of {n_gen} genomes published. Rows are in tree order; a heavy line starts each superfamily.
    Copies = firm + soft assigned; sim = median copy bitscore / consensus self-bitscore; green = 1,000 copies or more.
    MEG lists the four families (RL, RS, T2, TR). The tree badge is their sum.
    Prior evidence is the SINE label on the annotated copy of the Hao et al. figure.
    De novo = AnnoSINE2 + SINEbase-fragment scan candidates, clustered by SubFam for manual review (hla and tbr first).
  </p>
  <a class="btn secondary" href="chiroptera/LOG.md">Analysis Log</a>
  <a class="btn secondary" href="rhin.html">Rhin-1 peel alignments (<i>R. sinicus</i>)</a>
</section>

</main>
<script>
const RAW = 'https://raw.githubusercontent.com/Toki-bio/Tal/main/';
const MSA = 'https://toki-bio.github.io/MSA-viewer/';
function openMSA(relPath, title) {{
  const url = RAW + relPath;
  window.open(MSA + '?url=' + encodeURIComponent(url) + '&title=' + encodeURIComponent(title), '_blank');
}}
</script>
</body></html>
"""
    open(PAGE, "w", encoding="utf8", newline="\n").write(page)
    print(f"{n_done}/{n_gen} published")


if __name__ == "__main__":
    sys.exit(main())
