#!/usr/bin/env python3
"""tiePM review page for run13 (rsi-rfe): chiroptera/orth_rsi-rfe/tiepm_review/.

Inputs (copies of therioserver:/home/toki/sol_bank/rsi-rfe/run13_addfrag_e3e3273/work/ and run12_tiepm_ccdfeef/work/)
in one folder given as argv[1], with sub-folders r13/ (aln_rsi-rfe_rejected.aln.gz, aln_rsi-rfe_SINE.aln.gz,
stat_doubles_rsi-rfe, stat_multi_rsi-rfe, stat_rescue_rsi-rfe) and r12/ (orth_rsi-rfe.tsv, aln_rsi-rfe_rejected.aln.gz).
Writes data/aln_rsi-rfe_tiePM.aln.gz (306 tiePM + 2 SINE comparison alignments), data/tiepm.tsv and index.html.
Stat fields are ComPair.sh's (code e3e3273): row 1/2 = first/second genomic row of the alignment;
OneSINE/TwoSINE = esl-alipid % identity of row 1/2 to the reference in the SINE columns, OneTwo = row 1 vs row 2 there,
SO/ST = identical nt of row 1/2 with the reference, R1/R2 = residues of row 1/2 in the SINE columns,
LF = left-flank columns, RFlength = right-flank columns, LS1/LS2 = residues of row 1/2 left of the SINE start,
REFID = % of the reference's residues matched by the best genomic row (guard, not used for the verdict).
"""
import csv, gzip, html, io, os, re, sys

SRC = sys.argv[1]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "orth_rsi-rfe", "tiepm_review")
os.makedirs(os.path.join(OUT, "data"), exist_ok=True)
BUNDLE = "data/aln_rsi-rfe_tiePM.aln.gz"
LOC = re.compile(r'^(\S+?):(\d+)-(\d+)\(([+-])\)')


def species(chrom):
    return "rsi" if chrom.startswith("NC_142") else "rfe"


def read_bundle(path, keep):
    """{file name: list of text lines} for ##FILE blocks whose name passes keep()."""
    out, cur = {}, None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("##FILE "):
                name = line[7:].strip()
                cur = name if keep(name) else None
                if cur:
                    out[cur] = []
            elif line.startswith("##"):
                continue
            elif cur:
                out[cur].append(line)
    return out


def parse_fasta(lines):
    names, seqs = [], []
    for l in lines:
        if l.startswith(">"):
            names.append(l[1:]); seqs.append([])
        elif seqs:
            seqs[-1].append(l)
    return names, ["".join(s) for s in seqs]


def read_stats(paths, want):
    st = {}
    for p in paths:
        for line in open(p, encoding="utf8"):
            f = line.split()
            if len(f) < 3 or f[1] not in want:
                continue
            d = dict(x.split("=", 1) for x in f[2:] if "=" in x)
            st[(f[0].lstrip("./"), f[1])] = d
    return st


def iv(s):
    m = LOC.match(s)
    return (m.group(1), int(m.group(2)), int(m.group(3))) if m else None


def overlaps(a, b):
    return a and b and a[0] == b[0] and a[1] < b[2] and b[1] < a[2]


# ---- run13 tiePM alignments + stats
r13 = os.path.join(SRC, "r13")
tie = read_bundle(os.path.join(r13, "aln_rsi-rfe_rejected.aln.gz"), lambda n: n.endswith(".tiePM"))
stats = read_stats([os.path.join(r13, f) for f in ("stat_doubles_rsi-rfe", "stat_multi_rsi-rfe", "stat_rescue_rsi-rfe")],
                   {"tiePM", "SINE"})

# ---- run12 calls by locus overlap (orth table, then rejected bundle)
r12 = os.path.join(SRC, "r12")
r12calls = []  # (rsi iv, rfe iv, call)
for r in csv.DictReader(open(os.path.join(r12, "orth_rsi-rfe.tsv"), encoding="utf8"), delimiter="\t"):
    a = {r["species1"]: iv(r["locus1"]), r["species2"]: iv(r["locus2"])}
    r12calls.append((a.get("rsi"), a.get("rfe"), r["status"]))
for name, lines in read_bundle(os.path.join(r12, "aln_rsi-rfe_rejected.aln.gz"), lambda n: True).items():
    names, _ = parse_fasta(lines)
    a = {}
    for n in names[:-1]:
        v = iv(n)
        if v:
            a.setdefault(species(v[0]), v)
    r12calls.append((a.get("rsi"), a.get("rfe"), "rejected " + name.rsplit(".", 1)[-1]))
by_chrom = {}
for c in r12calls:
    if c[0]:
        by_chrom.setdefault(c[0][0], []).append(c)


def run12_call(rsi, rfe):
    hits = sorted({c[2] for c in by_chrom.get(rsi[0], []) if overlaps(c[0], rsi) and overlaps(c[1], rfe)})
    acc = [h for h in hits if not h.startswith("rejected")]
    return "/".join(acc or hits) or "none"


rows = []
for name, lines in tie.items():
    names, seqs = parse_fasta(lines)
    key = name[:-len(".tiePM")]
    s = stats[(key, "tiePM")]
    g1, g2, ref = names[0], names[1], names[-1]
    i1, i2 = iv(g1), iv(g2)
    sp1, sp2 = species(i1[0]), species(i2[0])
    LF, RFl = int(s["LF"]), int(s["RFlength"])
    L = len(seqs[0]); RF = L - RFl
    res = lambda q, a, b: len(q[a:b].replace("-", ""))
    RR = res(seqs[-1], LF, RF)
    R1, R2 = int(s["R1"]), int(s["R2"])
    assert (R1, R2) == (res(seqs[0], LF, RF), res(seqs[1], LF, RF)), name
    rf1, rf2 = res(seqs[0], RF, L), res(seqs[1], RF, L)
    d = dict(alignment=name, cluster=key.split(".")[0], nrows=len(names), reference=ref.split()[-1],
             RR=RR, LF=LF, RFlength=RFl, REFID=int(s["REFID"]), OneTwo=int(s["OneTwo"]), SL=int(s["SL"]),
             FL=int(s["FL"]), FR=int(s["FR"]))
    per = {sp1: dict(loc=g1.split("::")[0], R=R1, SO=int(s["SO"]), ID=int(s["OneSINE"]), LS=int(s["LS1"]), RS=rf1),
           sp2: dict(loc=g2.split("::")[0], R=R2, SO=int(s["ST"]), ID=int(s["TwoSINE"]), LS=int(s["LS2"]), RS=rf2)}
    assert set(per) == {"rsi", "rfe"}, name
    for sp in ("rsi", "rfe"):
        for k, v in per[sp].items():
            d[f"{sp}_{k}"] = v
    d["row1"] = sp1
    d["both_carry"] = int(R1 * 100 >= 50 * RR and R2 * 100 >= 50 * RR)
    d["same_LS"] = int(d["rsi_LS"] == d["rfe_LS"])
    d["run12"] = run12_call(iv(d["rsi_loc"]), iv(d["rfe_loc"]))
    rows.append(d)
rows.sort(key=lambda d: (-d["REFID"], d["alignment"]))
assert len(rows) == 306, len(rows)

# ---- comparison: 2 ordinary SINE alignments from the main bundle (3 rows, high REFID, OneTwo ~ 85-95, same junction)
sine_st = {k[0]: v for k, v in stats.items() if k[1] == "SINE"}
refs_hi = [d["reference"] for d in rows[:40]]
want_refs = []
for r in refs_hi:
    if r not in want_refs:
        want_refs.append(r)
cands = {}
main = read_bundle(os.path.join(r13, "aln_rsi-rfe_SINE.aln.gz"),
                   lambda n: n.endswith(".cl.dbl.SINE") and n[:-5] in sine_st
                   and 85 <= int(sine_st[n[:-5]]["OneTwo"]) <= 93 and int(sine_st[n[:-5]]["REFID"]) >= 85
                   and sine_st[n[:-5]]["LS1"] == sine_st[n[:-5]]["LS2"] and int(sine_st[n[:-5]]["FR"]) >= 85)
for n, lines in sorted(main.items()):
    names, _ = parse_fasta(lines)
    if len(names) == 3:
        cands.setdefault(names[-1].split()[-1], []).append(n)
comparison = []
for r in want_refs:
    if r in cands and len(comparison) < 2:
        comparison.append((cands[r][0], r))
assert len(comparison) == 2, (want_refs, list(cands)[:10])
comp_rows = []
for n, r in comparison:
    s = sine_st[n[:-5]]
    names, seqs = parse_fasta(main[n])
    sp1 = species(iv(names[0])[0])
    o1, o2 = (int(s["OneSINE"]), int(s["TwoSINE"]))
    ids = {sp1: o1, ("rfe" if sp1 == "rsi" else "rsi"): o2}
    comp_rows.append(dict(alignment=n, reference=r, REFID=int(s["REFID"]), OneTwo=int(s["OneTwo"]),
                          rsi_ID=ids["rsi"], rfe_ID=ids["rfe"], LS1=s["LS1"], LS2=s["LS2"]))

# ---- write bundle (same format as data/questioned/)
with gzip.open(os.path.join(OUT, BUNDLE), "wt", newline="\n") as fh:
    fh.write("##SINE_orth_loc bundle v1\n")
    for d in rows:
        fh.write(f"##FILE {d['alignment']}\n" + "\n".join(tie[d['alignment']]) + "\n")
    for n, _ in comparison:
        fh.write(f"##FILE {n}\n" + "\n".join(main[n]) + "\n")

cols = ["alignment", "cluster", "reference", "rsi_loc", "rfe_loc", "REFID", "rsi_R", "rfe_R", "RR", "rsi_SO", "rfe_SO",
        "rsi_ID", "rfe_ID", "OneTwo", "SL", "LF", "RFlength", "rsi_LS", "rfe_LS", "rsi_RS", "rfe_RS", "FL", "FR",
        "row1", "both_carry", "same_LS", "run12"]
with open(os.path.join(OUT, "data", "tiepm.tsv"), "w", encoding="utf8", newline="\n") as fh:
    fh.write("\t".join(cols) + "\n")
    for d in rows:
        fh.write("\t".join(str(d[c]) for c in cols) + "\n")

# ---- counts
def cnt(f):
    return sum(1 for d in rows if f(d))
bins = [("REFID < 50", lambda d: d["REFID"] < 50), ("50-59", lambda d: 50 <= d["REFID"] < 60),
        ("60-69", lambda d: 60 <= d["REFID"] < 70), ("&ge; 70", lambda d: d["REFID"] >= 70)]
bin_n = [(b, cnt(f)) for b, f in bins]
ref_n = {}
for d in rows:
    ref_n[d["reference"]] = ref_n.get(d["reference"], 0) + 1
r12_n = {}
for d in rows:
    r12_n[d["run12"]] = r12_n.get(d["run12"], 0) + 1
n_same = cnt(lambda d: d["same_LS"])
n_same_both = cnt(lambda d: d["same_LS"] and d["rsi_RS"] == d["rfe_RS"])
n_95 = cnt(lambda d: d["OneTwo"] >= 95)
n_95_same = cnt(lambda d: d["OneTwo"] >= 95 and d["same_LS"])
n_both = cnt(lambda d: d["both_carry"])
n_three = cnt(lambda d: d["nrows"] == 3)

# ---- look first: 4 high, 4 mid, 2 low (distinct references where possible), + 2 comparison SINE
def pick(pool, k):
    out, seen = [], set()
    for d in pool:
        if d["reference"] not in seen and d["nrows"] == 3:
            out.append(d); seen.add(d["reference"])
        if len(out) == k:
            return out
    for d in pool:
        if d not in out:
            out.append(d)
        if len(out) == k:
            break
    return out
hi = pick([d for d in rows if d["REFID"] >= 85 and d["same_LS"]], 4)
midpool = sorted([d for d in rows if 60 <= d["REFID"] < 70], key=lambda d: (abs(d["REFID"] - 65), d["alignment"]))
mid = pick(midpool, 4)
low = pick(sorted([d for d in rows if d["REFID"] < 50], key=lambda d: -d["REFID"]), 2)

PAGE = "https://toki-bio.github.io/Tal/chiroptera/orth_rsi-rfe/tiepm_review/"
def link(n):
    return f"../sine_loci_browser.html?bundle=tiepm_review/{BUNDLE}&name={n}&title={n}"
# browser resolves bundle relative to its own location (orth_rsi-rfe/), so bundle path is tiepm_review/data/...

def lf_li(d, why):
    return (f'<li><a href="{html.escape(link(d["alignment"]))}" target="_blank">{html.escape(d["alignment"])}</a> '
            f'&mdash; {html.escape(d["reference"])}, REFID {d["REFID"]}, rsi {d["rsi_ID"]} % / rfe {d["rfe_ID"]} % to the reference, '
            f'rsi vs rfe {d["OneTwo"]} %, junction {d["rsi_LS"]}/{d["rfe_LS"]}{why}</li>')
look = []
look += [lf_li(d, "") for d in hi]
look_mid = [lf_li(d, "") for d in mid]
look_low = [lf_li(d, "") for d in low]
look_cmp = [(f'<li><a href="{html.escape(link(c["alignment"]))}" target="_blank">{html.escape(c["alignment"])}</a> '
             f'&mdash; ordinary SINE call, {html.escape(c["reference"])}, REFID {c["REFID"]}, rsi {c["rsi_ID"]} % / rfe {c["rfe_ID"]} % '
             f'to the reference, rsi vs rfe {c["OneTwo"]} %, junction {c["LS1"]}/{c["LS2"]}</li>') for c in comp_rows]
look_urls = [PAGE.rsplit("tiepm_review/", 1)[0] + link(d["alignment"])[3:] for d in hi + mid + low] + \
            [PAGE.rsplit("tiepm_review/", 1)[0] + link(c["alignment"])[3:] for c in comp_rows]
with open(os.path.join(OUT, "data", "look_first_urls.txt"), "w", encoding="utf8", newline="\n") as fh:
    fh.write("\n".join(u.replace(" ", "%20") for u in look_urls) + "\n")


def tr(d):
    a = html.escape(link(d["alignment"]))
    return (f'<tr><td><a href="{a}" target="_blank">{html.escape(d["alignment"])}</a></td><td>{html.escape(d["reference"])}</td>'
            f'<td>{html.escape(d["rsi_loc"])}</td><td>{html.escape(d["rfe_loc"])}</td><td class="num"><b>{d["REFID"]}</b></td>'
            f'<td class="num">{d["rsi_R"]}</td><td class="num">{d["rfe_R"]}</td><td class="num">{d["RR"]}</td>'
            f'<td class="num">{d["rsi_SO"]}</td><td class="num">{d["rfe_SO"]}</td>'
            f'<td class="num">{d["rsi_ID"]}</td><td class="num">{d["rfe_ID"]}</td><td class="num">{d["OneTwo"]}</td>'
            f'<td class="num">{d["LF"]}</td><td class="num">{d["RFlength"]}</td>'
            f'<td class="num">{d["rsi_LS"]}/{d["rfe_LS"]}</td><td>{html.escape(d["run12"])}</td></tr>')


style = open(os.path.join(OUT, "..", "index.html"), encoding="utf8").read()
style = style[style.index("<style>"):style.index("</head>")]
nl = "\n"
base_n = {}
for k, v in ref_n.items():
    b = base_n.setdefault(k.split(":")[0], [0, 0])
    b[0] += v
    if ":" in k:
        b[1] += v
ref_txt = ", ".join(f"{html.escape(k)} {v[0]}" + (f" ({v[1]} of them a partial span)" if v[1] else "")
                    for k, v in sorted(base_n.items(), key=lambda kv: (-kv[1][0], kv[0])))
r12_txt = ", ".join(f"{html.escape(k)} {v}" for k, v in sorted(r12_n.items(), key=lambda kv: (-kv[1], kv[0])))
page = f"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>tiePM review rsi-rfe</title>
{style}</head>
<body>
<header>
  <h1>tiePM review: rsi-rfe, run13</h1>
  <div class="sub">306 alignments the pipeline set aside because both rows carry a copy &middot; run13_addfrag_e3e3273 &middot; <a href="../index.html" style="color:#fff;">back to the orthologous loci page</a></div>
</header>
<main>
<section class="card">
  <h2>What tiePM means, and the question</h2>
  <p style="margin:0;font-size:.92rem;">Each alignment has three rows: the rsi window, the rfe window and the reference consensus. Both flanks
  match, so the two windows are the same place in the two genomes. The pipeline first asks whether it is a SINE present in both
  (the two copies at least 65 % identical to each other and enough identical bases with the reference); these alignments failed that test.
  It then asks which genome has the insertion, by counting each row's residues in the SINE columns: a row carries the copy if it has
  &ge; 50 % of the reference's residues there. Here <b>both</b> the rsi and the rfe rows do ({n_both} of 306; in the other
  {306 - n_both} neither does), so the rule cannot say PM or MP, and the locus is rejected as <b>tiePM</b>.
  <b>The question for you: is this a SINE present in both genomes, only diverged?</b> If yes, these loci would become SINE
  (present in both) instead of being rejected.</p>
</section>

<section class="card">
  <h2>Look at these first</h2>
  <p style="margin:0 0 6px;font-size:.88rem;">Each link opens one alignment in the alignment browser. REFID = % of the reference's residues
  matched by the better genomic row; junction = residues of rsi/rfe left of the SINE start (equal = the copy starts at the same place in both windows).</p>
  <b style="font-size:.9rem;">High REFID (clear matches)</b><ol style="margin:4px 0 8px;font-size:.88rem;">{nl.join(look)}</ol>
  <b style="font-size:.9rem;">Middle REFID (60-69)</b><ol start="5" style="margin:4px 0 8px;font-size:.88rem;">{nl.join(look_mid)}</ol>
  <b style="font-size:.9rem;">Low REFID (&lt; 50)</b><ol start="9" style="margin:4px 0 8px;font-size:.88rem;">{nl.join(look_low)}</ol>
  <b style="font-size:.9rem;">For comparison: ordinary SINE calls (present in both) from the main bundle</b><ol start="11" style="margin:4px 0 0;font-size:.88rem;">{nl.join(look_cmp)}</ol>
</section>

<section class="card">
  <h2>Counts</h2>
  <div class="facts">
    <span><b>306</b> tiePM alignments ({len({d['cluster'] for d in rows})} clusters)</span>
    {"".join(f"<span>{b}: <b>{n}</b></span>" for b, n in bin_n)}
  </div>
  <p style="margin:4px 0;font-size:.88rem;"><b>By reference consensus:</b> {ref_txt}.</p>
  <p style="margin:4px 0;font-size:.88rem;"><b>Same junction</b> (rsi and rfe rows have the same number of residues left of the SINE start): <b>{n_same}</b> of 306
  ({n_same_both} of them also with the same number of residues right of the SINE).</p>
  <p style="margin:4px 0;font-size:.88rem;"><b>rsi and rfe rows &ge; 95 % identical to each other in the SINE columns</b> (shared-ancestral signal): <b>{n_95}</b> of 306
  ({n_95_same} of them with the same junction).</p>
  <p style="margin:4px 0;font-size:.88rem;"><b>run12 call of the same loci</b> (rsi and rfe windows overlapping a run12 row): {r12_txt}.</p>
</section>

<section class="card">
  <h2>All 306 tiePM alignments</h2>
  <p style="margin:0 0 8px;font-size:.84rem;" class="muted">Sorted by REFID, highest first. R = residues of the row in the SINE columns (ref = the reference's);
  SO = identical bases with the reference; id = % identity to the reference in the SINE columns (esl-alipid); rsi:rfe = % identity of the two rows there;
  LF/RF = left/right flank length in alignment columns; junction = residues left of the SINE start, rsi/rfe.</p>
  <p style="margin:10px 0;"><input id="qq" type="text" placeholder="filter: alignment, coordinate, reference or run12 call" style="width:340px;max-width:100%;padding:4px 6px;">
  <span id="qshown" style="color:var(--muted);font-size:.85rem;margin-left:8px;"></span>
  <a class="btn secondary" href="data/tiepm.tsv" style="margin-left:8px;">tiepm.tsv</a>
  <a class="btn secondary" href="{html.escape(f'../sine_loci_browser.html?bundle=tiepm_review/{BUNDLE}&title=rsi-rfe tiePM run13')}" target="_blank">open all 308</a></p>
  <div style="overflow:auto;max-height:720px;">
  <table class="tbl" id="loci">
    <thead><tr><th>Alignment (opens it)</th><th>Reference</th><th>rsi</th><th>rfe</th><th>REFID</th><th>R rsi</th><th>R rfe</th><th>R ref</th>
    <th>SO rsi</th><th>SO rfe</th><th>id rsi</th><th>id rfe</th><th>rsi:rfe</th><th>LF</th><th>RF</th><th>junction</th><th>run12 call</th></tr></thead>
    <tbody>
{nl.join(tr(d) for d in rows)}
    </tbody>
  </table>
  </div>
</section>
</main>
<script>
(function () {{
  const box = document.getElementById('qq'), shown = document.getElementById('qshown');
  const trs = Array.from(document.querySelectorAll('#loci tbody tr'));
  function apply() {{
    const t = box.value.trim().toLowerCase(); let n = 0;
    trs.forEach(tr => {{ const ok = !t || tr.textContent.toLowerCase().includes(t); tr.style.display = ok ? '' : 'none'; if (ok) n++; }});
    shown.textContent = n + ' of ' + trs.length + ' shown';
  }}
  box.addEventListener('input', apply); apply();
}})();
</script>
</body></html>
"""
open(os.path.join(OUT, "index.html"), "w", encoding="utf8", newline="\n").write(page)
print("bins", bin_n, "refs", ref_n, "run12", r12_n, "sameLS", n_same, n_same_both, "95", n_95, n_95_same,
      "both", n_both, "3rows", n_three, "comparison", [c["alignment"] for c in comp_rows])
print("\n".join(look_urls))
