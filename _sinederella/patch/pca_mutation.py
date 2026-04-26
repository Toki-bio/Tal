#!/usr/bin/env python3
"""
Mutation-space PCA for SINE copies.

Workflow:
  1. Stratified sample from assigned.fasta
  2. mafft --auto to align all consensus sequences
  3. mafft --add --keeplength to add sampled copies into consensus alignment
  4. Build binary mutation matrix: copy × variable_alignment_column
     (1 if copy nucleotide ≠ assigned-subfamily consensus, 0 otherwise)
  5. Filter to variable columns (mutation rate 3–97%)
  6. SVD → PC1/PC2 scatter coloured by assigned subfamily

Usage:
    python3 pca_mutation.py <run_root> [n_per_sf=200] [out.html=/tmp/pca_mutation.html]
"""

import subprocess, random, sys, json, shutil, tempfile
from pathlib import Path
from collections import defaultdict

try:
    import numpy as np
except ImportError:
    sys.exit("numpy required")

PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.2.min.js"
SF_PALETTE = [
    "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
    "#42d4f4", "#f032e6", "#bfef45", "#469990", "#9a6324",
    "#800000", "#aaffc3", "#808000", "#dcbeff", "#fabed4",
]


def read_fasta(path):
    out, name, seq = [], None, []
    with open(path, errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(seq)))
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())
    if name is not None:
        out.append((name, "".join(seq)))
    return out


def write_fasta(records, path):
    with open(path, "w") as fh:
        for n, s in records:
            fh.write(f">{n}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i+60] + "\n")


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: pca_mutation.py <run_root> [n_per_sf=200] [out.html]")

    run_root = Path(sys.argv[1])
    n_per_sf = int(sys.argv[2]) if len(sys.argv) > 2 else 200
    out_html = Path(sys.argv[3]) if len(sys.argv) > 3 else Path("/tmp/pca_mutation.html")

    # ── find step2 output ──────────────────────────────────────────────────
    s2 = None
    for cand in [run_root / "step2/step2_output",
                 run_root / "results/step2_output"]:
        if cand.is_dir():
            s2 = cand
            break
    if s2 is None:
        sys.exit(f"No step2 output under {run_root}")

    cons_fa = None
    for cand in [run_root / "consensuses.clean.fa",
                 run_root / "results/consensuses.fa"]:
        if cand.is_file():
            cons_fa = cand
            break
    if cons_fa is None:
        sys.exit("No consensuses FASTA found")

    assigned_fa = s2 / "assigned.fasta"
    if not assigned_fa.is_file():
        sys.exit(f"No assigned.fasta at {assigned_fa}")

    # ── read consensuses ───────────────────────────────────────────────────
    cons_records = read_fasta(cons_fa)
    cons_names = set(n for n, _ in cons_records)
    print(f"Consensuses ({len(cons_names)}): {sorted(cons_names)}")

    # ── read assigned.fasta, parse sf from header >id|sf|score ────────────
    print(f"Reading {assigned_fa} ...")
    all_copies = read_fasta(assigned_fa)
    seq_to_sf = {}
    for n, s in all_copies:
        parts = n.split("|")
        seq_to_sf[n] = parts[1] if len(parts) >= 2 else "unknown"

    # ── stratified sample ──────────────────────────────────────────────────
    random.seed(42)
    by_sf = defaultdict(list)
    for rec in all_copies:
        by_sf[seq_to_sf[rec[0]]].append(rec)

    sampled = []
    sf_order = sorted(by_sf.keys())
    for sf in sf_order:
        recs = by_sf[sf]
        take = min(len(recs), n_per_sf)
        sampled.extend(random.sample(recs, take))

    print(f"Sampled {len(sampled)} copies:")
    for sf in sf_order:
        cnt = sum(1 for n, _ in sampled if seq_to_sf[n] == sf)
        print(f"  {sf}: {cnt} of {len(by_sf[sf])}")

    # ── simplify copy names for mafft (avoid header parsing issues) ────────
    simple_to_sf = {}   # cp{i} → sf_name
    simple_records = []
    for i, (orig, seq) in enumerate(sampled):
        sname = f"cp{i}"
        simple_to_sf[sname] = seq_to_sf[orig]
        simple_records.append((sname, seq))

    tmpdir = Path(tempfile.mkdtemp(prefix="pca_mut_"))
    try:
        copies_fa = tmpdir / "copies.fa"
        cons_tmp  = tmpdir / "cons.fa"
        cons_aln  = tmpdir / "cons_aln.fa"
        all_aln   = tmpdir / "all_aln.fa"

        write_fasta(simple_records, copies_fa)
        write_fasta(cons_records,   cons_tmp)

        # ── Step 1: align consensus sequences ─────────────────────────────
        print(f"\nAligning {len(cons_records)} consensuses with mafft ...")
        r = subprocess.run(
            ["mafft", "--auto", "--quiet", "--thread", "8", "--nuc",
             str(cons_tmp)],
            capture_output=True, timeout=120)
        if r.returncode != 0:
            sys.exit(f"mafft (consensus) failed:\n{r.stderr.decode()[:500]}")
        cons_aln.write_bytes(r.stdout)
        cons_aln_recs = read_fasta(cons_aln)
        aln_len_cons = len(cons_aln_recs[0][1]) if cons_aln_recs else 0
        print(f"  Consensus alignment: {len(cons_aln_recs)} seqs × {aln_len_cons} cols")

        # ── Step 2: add copies (fixed length) ─────────────────────────────
        n_copies_total = len(simple_records)
        print(f"Adding {n_copies_total} copies (mafft --add --keeplength) ...")
        r2 = subprocess.run(
            ["mafft", "--add", str(copies_fa),
             "--keeplength", "--quiet", "--thread", "8", "--nuc",
             str(cons_aln)],
            capture_output=True, timeout=900)
        if r2.returncode != 0:
            sys.exit(f"mafft --add failed:\n{r2.stderr.decode()[:500]}")
        all_aln.write_bytes(r2.stdout)

        all_records = read_fasta(all_aln)
        print(f"  Full alignment: {len(all_records)} seqs")

        # ── parse alignment ────────────────────────────────────────────────
        cons_aln_dict = {}
        copy_aln      = []
        for name, seq in all_records:
            if name in cons_names:
                cons_aln_dict[name] = seq.upper()
            else:
                copy_aln.append((name, seq.upper()))

        print(f"  Parsed: {len(cons_aln_dict)} consensuses, {len(copy_aln)} copies")
        if not cons_aln_dict or not copy_aln:
            sys.exit("Alignment parsing failed")

        aln_len   = len(next(iter(cons_aln_dict.values())))
        n_copies  = len(copy_aln)
        copy_names = [n for n, _ in copy_aln]

        # ── build binary mutation matrix ───────────────────────────────────
        # 1 where copy ≠ assigned-subfamily consensus, 0 where they match or gap
        print(f"Building mutation matrix ({n_copies} × {aln_len}) ...")
        X_raw = np.zeros((n_copies, aln_len), dtype=np.float32)
        assigned_sfs = []

        for i, (cname, cseq) in enumerate(copy_aln):
            sf = simple_to_sf.get(cname, "unknown")
            assigned_sfs.append(sf)
            ref_seq = cons_aln_dict.get(sf) or next(iter(cons_aln_dict.values()))
            for j in range(aln_len):
                cc = cseq[j]   if j < len(cseq)    else '-'
                rc = ref_seq[j] if j < len(ref_seq) else '-'
                if cc not in '-Nn' and rc not in '-Nn' and cc != rc:
                    X_raw[i, j] = 1.0

        # filter variable columns: mutation rate 3 %–97 %
        col_rate   = X_raw.mean(axis=0)
        keep_mask  = (col_rate >= 0.03) & (col_rate <= 0.97)
        X          = X_raw[:, keep_mask]
        n_var      = int(keep_mask.sum())
        print(f"  Variable columns: {n_var} (of {aln_len})")

        if n_var < 2:
            sys.exit("Fewer than 2 variable columns — cannot run PCA")

        # per-subfamily average mutation rate at variable positions (diagnostic)
        print("\nPer-subfamily mutation profile (mean rate at variable cols):")
        for sf in sorted(set(assigned_sfs)):
            idx_sf = [i for i, s in enumerate(assigned_sfs) if s == sf]
            if idx_sf:
                mu = float(X[np.array(idx_sf), :].mean())
                print(f"  {sf}: {mu:.3f}  (n={len(idx_sf)})")

        # ── SVD PCA ───────────────────────────────────────────────────────
        print("\nRunning SVD ...")
        X_c = (X - X.mean(axis=0)).astype(np.float64)
        U, S_sv, Vt = np.linalg.svd(X_c, full_matrices=False)
        pc1       = (X_c @ Vt[0]).tolist()
        pc2       = (X_c @ Vt[1]).tolist()
        var_total = float((S_sv**2).sum()) or 1.0
        pct1      = round(float(S_sv[0]**2) / var_total * 100, 1)
        pct2      = round(float(S_sv[1]**2) / var_total * 100, 1)
        print(f"  PC1={pct1}%, PC2={pct2}%")

        # how much variance is captured by SF label (eta²)?
        sf_labels_num = np.array([sorted(set(assigned_sfs)).index(s)
                                   for s in assigned_sfs])
        pc1_arr = np.array(pc1)
        grand_mean = pc1_arr.mean()
        ss_total   = ((pc1_arr - grand_mean)**2).sum()
        ss_between = sum(
            (pc1_arr[sf_labels_num == k] - grand_mean).sum() ** 2 /
            (sf_labels_num == k).sum()
            for k in range(len(set(assigned_sfs)))
            if (sf_labels_num == k).sum() > 0
        )
        eta2 = ss_between / ss_total if ss_total > 0 else 0.0
        print(f"  eta² (subfamily → PC1): {eta2:.3f}  "
              f"(0=no separation, 1=perfect separation)")

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # ── build Plotly figure ────────────────────────────────────────────────
    sf_list   = sorted(set(assigned_sfs))
    traces    = []
    for idx, sf in enumerate(sf_list):
        xi, yi, ti = [], [], []
        for i, cname in enumerate(copy_names):
            if assigned_sfs[i] == sf:
                xi.append(round(float(pc1[i]), 4))
                yi.append(round(float(pc2[i]), 4))
                ti.append(cname)
        if not xi:
            continue
        traces.append({
            "type": "scatter", "mode": "markers",
            "x": xi, "y": yi, "name": sf, "text": ti,
            "marker": {
                "size": 5, "opacity": 0.7,
                "color": SF_PALETTE[idx % len(SF_PALETTE)],
            },
            "hovertemplate": (
                "<b>%{fullData.name}</b><br>%{text}<br>"
                "PC1=%{x:.3f}  PC2=%{y:.3f}<extra></extra>"),
        })

    layout = {
        "title": (
            f"Mutation-space PCA — mafft-aligned copies vs assigned-subfamily consensus"
            f"<br><sub>n={n_copies:,} copies &middot; {n_var} variable cols &middot; "
            f"PC1={pct1}% &middot; PC2={pct2}% variance &middot; "
            f"eta²(SF&#8594;PC1)={eta2:.3f}</sub>"
        ),
        "xaxis": {"title": f"PC1 ({pct1}% variance)", "zeroline": True},
        "yaxis": {"title": f"PC2 ({pct2}% variance)", "zeroline": True},
        "legend": {"title": {"text": "Assigned subfamily"}},
        "height": 720,
        "margin": {"t": 90, "r": 20, "b": 60, "l": 70},
    }

    fig_init = (
        f"Plotly.newPlot('plot', {json.dumps(traces)}, "
        f"{json.dumps(layout)}, {{responsive:true,displaylogo:false}});"
    )

    run_name = run_root.name
    html = f"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Mutation PCA &mdash; {run_name}</title>
<script src="{PLOTLY_CDN}" charset="utf-8"></script>
<style>
  body{{font-family:system-ui,sans-serif;margin:0;padding:20px;background:#f8f9fa}}
  h1{{font-size:1.25rem;margin:0 0 6px}}
  .sub{{font-size:.85rem;color:#555;margin-bottom:16px;line-height:1.5}}
  .key{{background:#f0f4ff;border-left:3px solid #4363d8;
        padding:8px 10px;font-size:.82rem;margin-bottom:14px;border-radius:2px}}
  .panel{{background:#fff;border:1px solid #ddd;border-radius:8px;padding:14px}}
</style>
</head><body>
<h1>Mutation-space PCA &mdash; SAQ Tal SINE</h1>
<div class="sub">
  Run: <code>{run_name}</code> &middot;
  {n_copies:,} copies ({n_per_sf} per subfamily) &middot;
  {len(sf_list)} subfamilies &middot;
  {n_var} variable alignment positions
</div>
<div class="key">
  <b>Method:</b> All sampled copies aligned to consensus bank via
  <code>mafft --add --keeplength</code>. Binary feature vector per copy:
  1 at each alignment column where the copy nucleotide &ne; the
  assigned-subfamily consensus (ignoring gaps and N).
  Variable columns only (mutation rate 3&ndash;97&thinsp;%).
  SVD on centred binary matrix.<br>
  <b>eta&sup2;(subfamily &rarr; PC1) = {eta2:.3f}</b>
  &nbsp;(0 = no separation, 1 = perfect; anything &gt; 0.1 is meaningful).
</div>
<div class="panel">
  <div id="plot"></div>
</div>
<script>{fig_init}</script>
</body></html>"""

    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html, encoding="utf-8")
    print(f"\nOutput: {out_html}  ({out_html.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
