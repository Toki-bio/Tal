#!/usr/bin/env python3
"""Add the SINE Tribes section to Chiroptera species pages (<code>/report.html).

Same section as the saq/ccr pages (_sinederella/patch/inject_tribes_section.py): metrics,
the 30 largest tribes, MSA links for the 10 aligned ones, the full TSV. The columns holding
hand annotations on ccr are left out; there are none yet for the bats.

Input per genome, made on therioserver by ~/chiro/chiro_tribes.sh (the eri extract_tribes.sh
method: vsearch --cluster_fast --id 0.99 --strand plus on the Step1 extracted pool, tribes of
>= 5 copies, top 10 aligned from up to 500 copies with 50 bp / 70 bp genomic flanks):
    <code>/tribes/<code>_tribes_summary.tsv
    <code>/tribes/<code>_tribeNN_<n>seqs.aln.fa
    <code>/tribes/tribes_metadata.txt

Idempotent: an existing section (between the SINE_TRIBES markers) is replaced.
Usage: inject_tribes.py [code ...]   (default: every genome in genomes.tsv with a tribes/ folder)
"""
import csv
import html
import sys
from pathlib import Path
from urllib.parse import quote

ROOT = Path(__file__).resolve().parents[1]
REPO_RAW_BASE = "https://raw.githubusercontent.com/Toki-bio/Tal/main"
VIEWER_BASE = "https://toki-bio.github.io/MSA-viewer/"
MAX_TABLE_ROWS = 30
MAX_ALIGNMENT_ROWS = 10
START = "  <!-- SINE_TRIBES_SECTION_START -->"
END = "  <!-- SINE_TRIBES_SECTION_END -->"


def fmt_int(v):
    return f"{int(v):,}"


def msa_link(code, filename, title, label):
    raw = f"{REPO_RAW_BASE}/{code}/tribes/{filename}"
    href = f"{VIEWER_BASE}?url={quote(raw, safe='')}&title={quote(title, safe='')}"
    return f"<a class='aln-link' href='{href}' target='_blank'>{html.escape(label)}</a>"


def alignment_file(tdir, code, tribe_id):
    m = sorted(tdir.glob(f"{code}_{tribe_id}_*seqs.aln.fa"))
    return m[0].name if m else None


def dominant(r):
    """Firm dominant subfamily; for a tribe with no firmly assigned member, the soft call.

    Soft = the copy failed 10/10 unanimity (or the bitscore bar) and was labelled by the query
    that found it (results/unassigned.tsv). The MEG-RS / MEG-RL overlap makes many such tribes."""
    if r["dominant_subfamily"] != "NA":
        return r["dominant_subfamily"]
    if r.get("soft_dominant"):
        return r["soft_dominant"] + " (soft)"
    return "unassigned"


def build_section(code, rows, tdir):
    total = len(rows)
    in_tribes = sum(int(r["seq_count"]) for r in rows)
    top_n = int(rows[0]["seq_count"]) if rows else 0
    top_sf = dominant(rows[0]) if rows else "n/a"
    tsv = f"{code}_tribes_summary.tsv"

    trs = []
    for r in rows[:MAX_TABLE_ROWS]:
        rank, tid = int(r["rank"]), r["tribe_id"]
        fn = alignment_file(tdir, code, tid) if rank <= MAX_ALIGNMENT_ROWS else None
        if fn:
            n_aln = fn.rsplit("_", 1)[1].replace("seqs.aln.fa", "")
            link = msa_link(code, fn, f"{code} {tid} tribe alignment", f"open ({n_aln} seqs)")
        else:
            link = "<span class='muted small'>not aligned</span>"
        dom = dominant(r)
        bd = r["subfamily_breakdown"].replace(",", "; ")
        if r.get("soft_breakdown"):
            bd = (bd + "; " if bd else "") + "soft " + r["soft_breakdown"].replace(",", "; soft ")
        bd = bd or "-"
        if len(bd) > 80:
            bd = bd[:77] + "..."
        trs.append(
            "<tr>"
            f"<td>{rank}</td>"
            f"<td><code>{html.escape(tid)}</code></td>"
            f"<td class='hl'>{fmt_int(r['seq_count'])}</td>"
            f"<td>{fmt_int(r['assigned_count'])}</td>"
            f"<td>{html.escape(dom)}</td>"
            f"<td>{fmt_int(r['dominant_subfamily_count'])}</td>"
            f"<td>{float(r['dominant_subfamily_pct_total']):.1f}%</td>"
            f"<td>{html.escape(bd)}</td>"
            f"<td>{link}</td>"
            "</tr>")
    body = "".join(trs) or "<tr><td colspan='9'>No tribes with at least 5 copies.</td></tr>"
    note = (f" Showing the largest {MAX_TABLE_ROWS} tribes here; the TSV contains all {fmt_int(total)}."
            if total > MAX_TABLE_ROWS else "")

    return (
        START + "\n"
        "  <section class='card' id='tribes'>"
        "<h2>SINE Tribes &mdash; 99% copy clusters</h2>"
        "<p class='intro'>SINEtribes-style clustering of the full Step1 extracted hit bank "
        "(<code>genome.clean_step1/extracted.fasta</code>) using <code>vsearch --cluster_fast --id 0.99 --strand plus</code>, "
        "independent of subfamily assignment; the dominant subfamily is read from the assignment afterwards "
        "(<i>soft</i>: members that failed the 10/10 vote, labelled by the query that found them). "
        "Tribes require at least 5 copies. The 10 largest are aligned from up to 500 copies (random sample, seed 42), "
        "each re-extracted with 50 bp upstream and 70 bp downstream genomic flanks from <code>genome.clean.fa</code> "
        "and aligned with MAFFT (<code>--localpair --maxiterate 1000</code>). "
        f"<a class='aln-link green' href='tribes/{html.escape(tsv)}' target='_blank'>Full TSV</a>{note}</p>"
        "<div class='metrics'>"
        f"<div class='metric'><div class='num'>{fmt_int(total)}</div><div class='lbl'>Tribes found</div></div>"
        f"<div class='metric'><div class='num'>{fmt_int(in_tribes)}</div><div class='lbl'>Copies in tribes</div></div>"
        f"<div class='metric'><div class='num'>{fmt_int(top_n)}</div><div class='lbl'>Largest tribe copies</div></div>"
        f"<div class='metric'><div class='num'>{html.escape(top_sf)}</div><div class='lbl'>Largest tribe subfamily</div></div>"
        "</div>"
        "<table class='tbl'><thead><tr>"
        "<th>Rank</th><th>Tribe</th><th class='hl'>Sequences</th><th>Assigned</th><th>Dominant subfamily</th>"
        "<th>Dominant count</th><th>Dominant % total</th><th>Assigned breakdown</th><th>Alignment</th>"
        "</tr></thead><tbody>"
        f"{body}"
        "</tbody></table>"
        "</section>\n"
        + END)


def strip_existing(text):
    s = text.find(START)
    if s == -1:
        return text
    e = text.find(END, s)
    if e == -1:
        raise ValueError("SINE_TRIBES_SECTION_START without END")
    e += len(END)
    while e < len(text) and text[e] == "\n" and e - (text.find(END, s) + len(END)) < 2:
        e += 1
    return text[:s] + text[e:]


def patch(code):
    rep = ROOT / code / "report.html"
    tdir = ROOT / code / "tribes"
    summ = tdir / f"{code}_tribes_summary.tsv"
    if not rep.exists() or not summ.exists():
        return f"{code}: skipped (no report or no tribes summary)"
    with summ.open(newline="", encoding="utf8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    text = strip_existing(rep.read_text(encoding="utf8"))
    if '<a href="#tribes">Tribes</a>' not in text:
        text = text.replace('    <a href="#alignments">Alignments</a>\n',
                            '    <a href="#alignments">Alignments</a>\n    <a href="#tribes">Tribes</a>\n', 1)
    marker = '  <section class="card" id="overview">'
    if marker not in text:
        raise ValueError(f"{rep}: no overview section to insert before")
    text = text.replace(marker, build_section(code, rows, tdir) + "\n\n" + marker, 1)
    rep.write_text(text, encoding="utf8", newline="\n")
    n_aln = len(list(tdir.glob(f"{code}_tribe*_*seqs.aln.fa")))
    return f"{code}: {len(rows)} tribes, {n_aln} alignments, largest {rows[0]['seq_count'] if rows else 0}"


def main(argv):
    codes = argv[1:]
    if not codes:
        with (ROOT / "chiroptera" / "genomes.tsv").open(encoding="utf8") as fh:
            codes = [r["code"] for r in csv.DictReader(fh, delimiter="\t")]
    for c in codes:
        print(patch(c))


if __name__ == "__main__":
    main(sys.argv)
