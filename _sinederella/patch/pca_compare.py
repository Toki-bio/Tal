#!/usr/bin/env python3
"""
Compare two PCA approaches for SINE copy visualization:
  A) SVD on column-normalised bitscore matrix (current step6_report.py)
  B) SINEplot approach: MDS of subfamily centres + weighted centroid placement

Usage:
    python3 pca_compare.py <run_root> <n_copies> <out_html>

Example:
    python3 pca_compare.py \
        /data/W/toki/Genomes/Mammalia/Eulipotyphla/saq/run_20260425_182219 \
        1500 /tmp/pca_compare.html
"""

import subprocess, random, tempfile, json, sys, shutil
from pathlib import Path
from collections import defaultdict

try:
    import numpy as np
except ImportError:
    sys.exit("numpy required")

PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.2.min.js"
SF_PALETTE = [
    "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
    "#42d4f4", "#f032e6", "#bfef45", "#fabed4", "#469990",
    "#dcbeff", "#9a6324", "#800000", "#aaffc3", "#808000",
]


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

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


def run_ssearch(query_fa, db_fa, threads=8, timeout=600):
    """Returns dict: (query, subject) -> max bitscore."""
    r = subprocess.run(
        ["ssearch36", "-m", "8", "-T", str(threads),
         str(query_fa), str(db_fa)],
        capture_output=True, timeout=timeout
    )
    scores = {}
    for line in r.stdout.decode(errors="replace").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        p = line.split("\t")
        if len(p) < 12:
            continue
        try:
            q, s, bs = p[0], p[1], float(p[11])
        except (ValueError, IndexError):
            continue
        key = (q, s)
        if key not in scores or scores[key] < bs:
            scores[key] = bs
    return scores


# ---------------------------------------------------------------------------
# Classical MDS
# ---------------------------------------------------------------------------

def classical_mds(D, n_components=2):
    """Classical (metric) MDS from square distance matrix D."""
    n = D.shape[0]
    D2 = D ** 2
    row_mean = D2.mean(axis=1, keepdims=True)
    col_mean = D2.mean(axis=0, keepdims=True)
    grand_mean = D2.mean()
    B = -0.5 * (D2 - row_mean - col_mean + grand_mean)
    eigvals, eigvecs = np.linalg.eigh(B)
    # sort descending
    idx = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[idx], eigvecs[:, idx]
    coords = eigvecs[:, :n_components] * np.sqrt(np.maximum(eigvals[:n_components], 0))
    return coords


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) < 4:
        sys.exit("usage: pca_compare.py <run_root> <n_copies> <out_html>")

    run_root = Path(sys.argv[1])
    n_copies_target = int(sys.argv[2])
    out_html = Path(sys.argv[3])

    s2 = None
    for cand in [run_root / "step2/step2_output",
                 run_root / "results/step2_output"]:
        if cand.is_dir():
            s2 = cand
            break
    if s2 is None:
        sys.exit(f"No step2 output found under {run_root}")

    cons_fa = None
    for cand in [run_root / "consensuses.clean.fa",
                 run_root / "results/consensuses.fa"]:
        if cand.is_file():
            cons_fa = cand
            break
    if cons_fa is None:
        sys.exit("No consensuses FASTA found")

    assigned_fa = s2 / "assigned.fasta"
    self_bits_tsv = s2 / "self_bits_real.tsv"

    # Read consensuses
    cons_records = read_fasta(cons_fa)
    cons_names = [n for n, _ in cons_records]
    n_cons = len(cons_names)
    print(f"Consensuses ({n_cons}): {cons_names}")

    # Read self-bits from file
    self_bits = {}
    if self_bits_tsv.is_file():
        with open(self_bits_tsv) as fh:
            for line in fh:
                p = line.strip().split("\t")
                if len(p) >= 2:
                    try:
                        self_bits[p[0]] = float(p[1])
                    except ValueError:
                        pass

    # Read assigned FASTA; parse subfamily from header ">id|sf|score"
    print(f"Reading {assigned_fa} ...")
    all_copies = read_fasta(assigned_fa)
    seq_to_sf = {}
    for n, s in all_copies:
        parts = n.split("|")
        seq_to_sf[n] = parts[1] if len(parts) >= 2 else "unknown"

    # Stratified sample
    random.seed(42)
    by_sf = defaultdict(list)
    for rec in all_copies:
        by_sf[seq_to_sf[rec[0]]].append(rec)

    sampled = []
    per_sf = max(1, n_copies_target // len(by_sf))
    for sf in sorted(by_sf):
        recs = by_sf[sf]
        take = min(len(recs), per_sf)
        sampled.extend(random.sample(recs, take))
    if len(sampled) > n_copies_target:
        sampled = random.sample(sampled, n_copies_target)

    print(f"Sampled {len(sampled)} copies across {len(by_sf)} subfamilies:")
    for sf in sorted(by_sf):
        cnt = sum(1 for n, _ in sampled if seq_to_sf[n] == sf)
        print(f"  {sf}: {cnt}")

    tmpdir = Path(tempfile.mkdtemp(prefix="pca_cmp_"))
    copies_tmp = tmpdir / "copies.fa"
    cons_tmp = tmpdir / "cons.fa"
    write_fasta(sampled, copies_tmp)
    write_fasta(cons_records, cons_tmp)

    copy_names = [n for n, _ in sampled]
    n_copies = len(copy_names)

    # ==========================================================================
    # SHARED: copies vs consensuses (both approaches use these scores)
    # ==========================================================================
    print(f"\nRunning ssearch36: {n_copies} copies × {n_cons} consensuses ...")
    scores_cv = run_ssearch(copies_tmp, cons_tmp)
    print(f"  → {len(scores_cv)} alignments")

    # ==========================================================================
    # APPROACH A: SVD PCA
    # ==========================================================================
    print("\n--- Approach A: SVD PCA ---")

    X = np.zeros((n_copies, n_cons), dtype=np.float32)
    for i, cname in enumerate(copy_names):
        for j, sname in enumerate(cons_names):
            sb = self_bits.get(sname, 1.0) or 1.0
            X[i, j] = scores_cv.get((cname, sname), 0.0) / sb

    X_c = (X - X.mean(axis=0)).astype(np.float64)
    U, S_sv, Vt = np.linalg.svd(X_c, full_matrices=False)
    pc1_A = (X_c @ Vt[0]).tolist()
    pc2_A = (X_c @ Vt[1]).tolist()
    var_total = float((S_sv**2).sum()) or 1.0
    pct1_A = round(float(S_sv[0]**2) / var_total * 100, 1)
    pct2_A = round(float(S_sv[1]**2) / var_total * 100, 1)
    print(f"  PC1={pct1_A}%, PC2={pct2_A}%")

    # ==========================================================================
    # APPROACH B: SINEplot MDS + weighted centroid
    # ==========================================================================
    print("\n--- Approach B: SINEplot MDS + weighted centroid ---")

    # Consensus vs consensus (for MDS distance matrix)
    print(f"  Running ssearch36: {n_cons} consensuses × {n_cons} consensuses ...")
    scores_cc = run_ssearch(cons_tmp, cons_tmp, timeout=60)
    print(f"  → {len(scores_cc)} alignments")

    # Self-scores: prefer ssearch36 self-alignment, fall back to file
    self_scores = {}
    for ci in cons_names:
        ss = scores_cc.get((ci, ci)) or self_bits.get(ci)
        self_scores[ci] = ss or 1.0

    # Distance matrix: 1 - cross_score / sqrt(self_i * self_j)
    D = np.zeros((n_cons, n_cons))
    for i, ci in enumerate(cons_names):
        for j, cj in enumerate(cons_names):
            if i == j:
                D[i, j] = 0.0
            else:
                cross = max(
                    scores_cc.get((ci, cj), 0.0),
                    scores_cc.get((cj, ci), 0.0),
                )
                norm = (self_scores[ci] * self_scores[cj]) ** 0.5 or 1.0
                D[i, j] = max(0.0, 1.0 - cross / norm)

    print("  Subfamily distance matrix (rounded):")
    for i, ci in enumerate(cons_names):
        row = "  " + ci.ljust(16) + " ".join(f"{D[i,j]:.3f}" for j in range(n_cons))
        print(row)

    # Classical MDS → 2D positions for subfamily centres
    sf_pos_2d = classical_mds(D, n_components=2)
    sf_positions = {cons_names[i]: sf_pos_2d[i].tolist() for i in range(n_cons)}
    print("  Subfamily MDS centres:")
    for sf, pos in sf_positions.items():
        print(f"    {sf}: ({pos[0]:.3f}, {pos[1]:.3f})")

    # Weighted centroid placement per copy
    x_B, y_B = [], []
    conflict_ratios = []
    for cname in copy_names:
        weights = {sf: max(0.0, scores_cv.get((cname, sf), 0.0))
                   for sf in cons_names}
        total_w = sum(weights.values())
        if total_w > 0:
            x = sum(weights[sf] * sf_positions[sf][0] for sf in cons_names) / total_w
            y = sum(weights[sf] * sf_positions[sf][1] for sf in cons_names) / total_w
        else:
            x = float(np.mean([p[0] for p in sf_positions.values()]))
            y = float(np.mean([p[1] for p in sf_positions.values()]))
        x_B.append(x)
        y_B.append(y)

        # conflict ratio = 2nd_best / 1st_best
        sorted_w = sorted(weights.values(), reverse=True)
        if len(sorted_w) >= 2 and sorted_w[0] > 0:
            conflict_ratios.append(sorted_w[1] / sorted_w[0])
        else:
            conflict_ratios.append(0.0)

    print(f"  Conflict ratio: mean={float(np.mean(conflict_ratios)):.3f}  "
          f"median={float(np.median(conflict_ratios)):.3f}")

    shutil.rmtree(tmpdir, ignore_errors=True)

    # ==========================================================================
    # Build HTML with both plots
    # ==========================================================================
    subfams_ordered = sorted(set(seq_to_sf[n] for n in copy_names))

    def make_traces(x_vals, y_vals):
        traces = []
        for idx, sf in enumerate(subfams_ordered):
            xi, yi, ti = [], [], []
            for i, cname in enumerate(copy_names):
                if seq_to_sf.get(cname) == sf:
                    xi.append(round(float(x_vals[i]), 4))
                    yi.append(round(float(y_vals[i]), 4))
                    ti.append(cname)
            if not xi:
                continue
            traces.append({
                "type": "scatter", "mode": "markers",
                "x": xi, "y": yi, "name": sf, "text": ti,
                "marker": {
                    "size": 5, "opacity": 0.65,
                    "color": SF_PALETTE[idx % len(SF_PALETTE)],
                },
                "hovertemplate": (
                    "<b>%{fullData.name}</b><br>%{text}<br>"
                    "x=%{x:.3f} y=%{y:.3f}<extra></extra>"),
            })
        return traces

    # Subfamily centre markers for plot B
    center_trace = {
        "type": "scatter", "mode": "markers+text",
        "x": [round(sf_positions[sf][0], 3) for sf in cons_names],
        "y": [round(sf_positions[sf][1], 3) for sf in cons_names],
        "text": cons_names,
        "textposition": "top center",
        "name": "Subfamily centres",
        "marker": {"size": 14, "symbol": "diamond", "color": "black", "opacity": 0.85},
        "showlegend": True,
        "hovertemplate": "<b>%{text}</b><br>MDS centre<extra></extra>",
    }

    fig_A = {
        "data": make_traces(pc1_A, pc2_A),
        "layout": {
            "title": (f"Approach A — SVD PCA on normalised bitscore matrix"
                      f"<br><sub>n={n_copies:,} copies &middot; "
                      f"PC1={pct1_A}% var &middot; PC2={pct2_A}% var</sub>"),
            "xaxis": {"title": f"PC1 ({pct1_A}% variance)", "zeroline": True},
            "yaxis": {"title": f"PC2 ({pct2_A}% variance)", "zeroline": True},
            "legend": {"title": {"text": "Subfamily"}},
            "height": 620,
            "margin": {"t": 80, "r": 20, "b": 60, "l": 70},
        },
    }

    fig_B = {
        "data": make_traces(x_B, y_B) + [center_trace],
        "layout": {
            "title": (f"Approach B — SINEplot: MDS centres + bitscore-weighted centroid"
                      f"<br><sub>n={n_copies:,} copies &middot; "
                      f"conflict ratio mean={float(np.mean(conflict_ratios)):.3f}</sub>"),
            "xaxis": {"title": "Dim 1 (MDS)", "zeroline": True},
            "yaxis": {
                "title": "Dim 2 (MDS)", "zeroline": True,
                "scaleanchor": "x", "scaleratio": 1,
            },
            "legend": {"title": {"text": "Subfamily"}},
            "height": 620,
            "margin": {"t": 80, "r": 20, "b": 60, "l": 70},
        },
    }

    fig_init = (
        f"Plotly.newPlot('plotA', {json.dumps(fig_A['data'])}, "
        f"{json.dumps(fig_A['layout'])}, {{responsive:true,displaylogo:false}});\n"
        f"Plotly.newPlot('plotB', {json.dumps(fig_B['data'])}, "
        f"{json.dumps(fig_B['layout'])}, {{responsive:true,displaylogo:false}});\n"
    )

    html = f"""<!doctype html>
<html lang="en"><head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>PCA comparison &mdash; SAQ SINEs</title>
<script src="{PLOTLY_CDN}" charset="utf-8"></script>
<style>
  body {{ font-family: system-ui, sans-serif; margin: 0; padding: 20px; background: #f8f9fa; }}
  h1 {{ font-size: 1.3rem; margin: 0 0 4px; }}
  .sub {{ font-size: .85rem; color: #666; margin-bottom: 16px; }}
  .grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }}
  @media (max-width: 1000px) {{ .grid {{ grid-template-columns: 1fr; }} }}
  .panel {{ background: #fff; border: 1px solid #ddd; border-radius: 8px; padding: 14px; }}
  .panel h2 {{ font-size: 1rem; margin: 0 0 6px; }}
  .desc {{ font-size: .82rem; color: #555; margin: 0 0 10px; line-height: 1.5; }}
  .key {{ background: #f0f4ff; border-left: 3px solid #4363d8;
          padding: 8px 10px; font-size: .82rem; margin-bottom: 10px; }}
</style>
</head><body>
<h1>SAQ SINE &mdash; PCA approach comparison</h1>
<div class="sub">
  {n_copies:,} copies stratified by subfamily &middot;
  {n_cons} consensuses &middot; ssearch36 -m 8 &middot; numpy SVD / MDS &middot;
  run: {run_root.name}
</div>
<div class="grid">
  <div class="panel">
    <h2>A &mdash; SVD PCA (current approach)</h2>
    <p class="desc">
      Each copy &rarr; vector of <code>bitscore / self-bits</code> vs each consensus
      &rarr; column-centred matrix &rarr; SVD &rarr; PC1, PC2.<br>
      Axes = directions of maximum <em>global</em> variance.
      No explicit subfamily positions; clusters appear only if variance is
      driven by subfamily identity.
    </p>
    <div id="plotA"></div>
  </div>
  <div class="panel">
    <h2>B &mdash; SINEplot: MDS + weighted centroid</h2>
    <p class="desc">
      Subfamily centres placed in 2D by <em>classical MDS</em> on a pairwise
      distance matrix (derived from inter-subfamily ssearch36 scores).
      Each copy placed at bitscore-weighted centroid of those centres.<br>
      Diamonds = subfamily centres. Copies are pulled toward their best match,
      so subfamily clouds form even when bitscores are similar.
    </p>
    <div id="plotB"></div>
  </div>
</div>
<script>{fig_init}</script>
</body></html>"""

    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html, encoding="utf-8")
    print(f"\nOutput: {out_html}  ({out_html.stat().st_size / 1024:.0f} KB)")


if __name__ == "__main__":
    main()
