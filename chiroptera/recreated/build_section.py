#!/usr/bin/env python3
"""Build the "Rebuilt consensi across species" card for chiroptera.html from the files in this folder
(recons.py / groups.py / nearest.py output, 2026-09-28) and insert it after the family-tree card.
Genome copy counts are read from the family tree's hover titles on the same page."""
import html, os, re
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
PAGE = os.path.join(HERE, "..", "..", "chiroptera.html")
MARK_A, MARK_B = "<!-- recreated-consensi:start -->", "<!-- recreated-consensi:end -->"


def tsv(name):
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(HERE, name), encoding="utf8") if l.strip()]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


page = open(PAGE, encoding="utf8").read()
# genome copy counts per (code, family) from the tree titles: "Rhinolophidae Rhin-1: rsi: 30,079 copies, sim 0.48; ..."
copies = {}
for fam_part, body in re.findall(r"<title>[A-Za-z]+ (Rhin-1|VES|MEG): ([^<]+)</title>", page):
    if fam_part == "MEG":
        for code, rest in re.findall(r"(\w{3}): ((?:(?:RL|RS|T2|TR): [\d,]+ copies, sim [\d.]+;? ?)+)", body):
            for sub, n in re.findall(r"(RL|RS|T2|TR): ([\d,]+) copies", rest):
                copies[(code, "MEG-" + sub)] = n
    else:
        for code, n in re.findall(r"(\w{3}): ([\d,]+) copies", body):
            copies[(code, fam_part)] = n

near = {r["name"]: r for r in tsv("nearest.tsv")}
groups = {r["name"]: r for r in tsv("groups.tsv")}
members = tsv("members.tsv")
skipped = tsv("skipped.tsv")

# groups at identity >= 0.80, with their members
g80 = OrderedDict()
for r in tsv("groups.tsv"):
    g = r["group_0.80"]
    if g != "-":
        g80.setdefault(g, []).append(r["name"])


def label(name):
    if name.startswith("REF_"):
        return "<b>%s</b> (query)" % html.escape(name[4:])
    m = near.get(name) or {}
    code, fam = m.get("code", "?"), m.get("sine_family", "?")
    c = copies.get((code, fam))
    return "%s %s%s" % (html.escape(code), html.escape(fam), (" <span style=\"color:var(--muted)\">(%s)</span>" % c) if c else "")


def gsummary(names):
    refs = sorted(n[4:] for n in names if n.startswith("REF_"))
    own = sorted({(near[n]["sine_family"]) for n in names if n in near})
    return refs, own


rows_g = []
for g, names in sorted(g80.items(), key=lambda x: -len(x[1])):
    refs, own = gsummary(names)
    what = ("contains the query " + ", ".join(refs)) if refs else "<b>no query in this group</b> - a different element"
    rows_g.append("<tr><td>%s</td><td class=\"num\">%d</td><td>%s</td><td>%s</td><td style=\"font-size:.82rem;\">%s</td></tr>"
                  % (g, len(names), ", ".join(html.escape(x) for x in own) or "-", what,
                     "; ".join(label(n) for n in sorted(names, key=lambda n: (n.startswith("REF_") is False, n)))))

rows_m = []
for m in sorted(members, key=lambda r: (r["sine_family"], -float(near.get(r["name"], {}).get("identity_to_own_query") or 0))):
    n = near.get(m["name"], {})
    own_id = n.get("identity_to_own_query", "")
    close = n.get("closest_ref", "-")
    odd = close != m["sine_family"] and not (m["sine_family"].startswith("r") and close.startswith("r"))
    style = " style=\"background:#fff4e5;\"" if odd else ""
    rows_m.append("<tr%s><td>%s</td><td>%s</td><td>%s</td><td class=\"num\">%s</td><td class=\"num\">%s</td><td class=\"num\">%s</td>"
                  "<td class=\"num\">%s</td><td>%s</td><td class=\"num\">%s</td><td>%s</td></tr>"
                  % (style, html.escape(m["sine_family"]), html.escape(m["code"]), html.escape(m["bat_family"]),
                     copies.get((m["code"], m["sine_family"]), ""), m["n"], "%.0f" % float(m["score"] or 0), m["len"],
                     html.escape(close) + (" <b>&#9888;</b>" if odd else ""), n.get("identity", ""),
                     own_id + ("" if groups.get(m["name"], {}).get("group_0.80", "-") == "-" else " &middot; " + groups[m["name"]]["group_0.80"])))

from collections import Counter
_sk = Counter(r["family"] for r in skipped)
sk = ", ".join("%s in %d species" % (f, k) for f, k in _sk.most_common())
n_core = sum(1 for m in members if m["name"].endswith("_core"))

card = """%s
<section class="card">
  <h2>Rebuilt consensi across species &mdash; what each query actually finds</h2>
  <p style="margin:0 0 10px;font-size:.9rem;color:var(--muted);">Row 1 of every top100 plate is the consensus <b>rebuilt from that species' own copies</b>.
  All %d of them (25 species, every SINE family with an element on its plate) are aligned here together with the %d queries as searched
  (<code>REF_</code>: SINEbase Rhin-1, VES, MEG-RL/RS/T2/TR and the ten <i>R. sinicus</i> subfamilies r1&ndash;r10).
  Taken from the plates re-run through the current chain (2026-09-28): the element as the verdict judges it &mdash; no trim proposals,
  additions only where the copies carry them; lowercase = a proposed addition past the original.
  Names: <code>family__code_batfamily_nCOPIES-ON-PLATE_sSCORE</code>; <code>_core</code> = an array plate (score &lt; 50, rebuilt row &gt; 2.5&times; the query)
  cut to the part over the query (%d rows); <code>_R_</code> = reverse-complemented by MAFFT (8 rows, all weak MEG plates).
  MAFFT <code>--localpair --maxiterate 1000 --adjustdirectionaccurately --reorder</code>, so similar consensi sit together.
  Not included (no element on the plate, or &lt; 5 copies): %d plates &mdash; %s.</p>
  <div class="actions" style="margin-bottom:10px;">
    <a class="btn" href="javascript:void(0)" onclick="openMSA('chiroptera/recreated/recreated.aln.fa','Rebuilt consensi, 25 bat species + queries')">Open the alignment (%d rows)</a>
    <a class="btn secondary" href="chiroptera/recreated/recreated.fa">unaligned FASTA</a>
    <a class="btn secondary" href="chiroptera/recreated/nearest.tsv">closest query per consensus</a>
    <a class="btn secondary" href="chiroptera/recreated/groups.tsv">groups</a>
    <a class="btn secondary" href="chiroptera/recreated/skipped.tsv">not included</a>
  </div>
  <h3 style="margin:14px 0 4px;font-size:1rem;">Groups &mdash; rebuilt consensi at &ge; 80 %% identity to one another (single linkage)</h3>
  <p style="margin:0 0 6px;font-size:.85rem;color:var(--muted);">Identity over the columns both have a letter (&ge; 60 of them and &ge; 60 %% of the shorter).
  Genome copy counts (from the tree above) in brackets.</p>
  <div style="overflow-x:auto;"><table class="tbl">
    <thead><tr><th>group</th><th>rows</th><th>searched as</th><th>query in group?</th><th>members</th></tr></thead>
    <tbody>
%s
    </tbody></table></div>
  <h3 style="margin:14px 0 4px;font-size:1rem;">Every rebuilt consensus &mdash; its closest query</h3>
  <p style="margin:0 0 6px;font-size:.85rem;color:var(--muted);">Shaded, &#9888;: the closest query is not the one it was found with.
  Identity to the closest query and to its own query, from the alignment; &middot; G = group above.</p>
  <div style="overflow-x:auto;"><table class="tbl">
    <thead><tr><th>searched as</th><th>species</th><th>bat family</th><th>genome copies</th><th>on plate</th><th>score</th><th>length</th>
    <th>closest query</th><th>identity</th><th>identity to own query &middot; group</th></tr></thead>
    <tbody>
%s
    </tbody></table></div>
</section>
%s""" % (MARK_A, len(members), sum(1 for l in open(os.path.join(HERE, "recreated.fa")) if l.startswith(">REF_")),
         n_core, len(skipped), html.escape(sk), len(members) + 16, "\n".join(rows_g), "\n".join(rows_m), MARK_B)

if MARK_A in page:
    page = re.sub(re.escape(MARK_A) + ".*?" + re.escape(MARK_B), lambda _: card, page, flags=re.S)
else:
    anchor = "</section>\n\n<section class=\"card\">\n  <h2>Pteropodoidea"
    assert page.count(anchor) == 1, "anchor after the family-tree card not found"
    page = page.replace(anchor, "</section>\n\n" + card + "\n\n<section class=\"card\">\n  <h2>Pteropodoidea", 1)
open(PAGE, "w", encoding="utf8", newline="\n").write(page)
print("card written: %d groups, %d rows, copies known for %d" % (len(rows_g), len(rows_m), sum(1 for m in members if (m["code"], m["sine_family"]) in copies)))
