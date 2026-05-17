"""
Build figures/conservation_figure.html — combined divergence-from-consensus
figure with 6 species panels (A–F), interactive reorder, caption, methods,
legend and per-panel PNG download.

Also writes figures/plot_data_{sp}.json for offline kaleido export.

Run from c:\work\Tal:
    .\.venv\Scripts\python.exe _make_figure_page.py
"""

import re, json, os, html as htmllib

# ---------------------------------------------------------------------------
# 1. Authoritative rename maps  (Table 1 v2.docx — verified)
# ---------------------------------------------------------------------------
RENAME = {
    'saq': {
        's1_30seqs':   'c_Saq (s1)',
        's2_38seqs':   'a_Saq',
        's3_43seqs':   'b1_Saq',
        's4_60seqs':   'd_Saq',
        's5_5seqs':    'e_Saq',
        's6_7seqs':    'f_Saq',
        's7g_172seqs': 'b2_Saq',
        's8_225seqs':  'g_Saq',
        's9_20seqs':   'c_Saq (s9)',
    },
    'ccr': {
        'g1_180seqs': 'b1_Ccr',
        'g2_103seqs': 'd_Ccr',
        'g3_71seqs':  'e_Ccr',
        'g4_77seqs':  'f_Ccr',
        'g5_7seqs':   'g_Ccr',
        'g6_58seqs':  'c_Ccr',
        'g7_3seqs':   'b2_Ccr',
    },
    'toc': {
        't1_45seqs':  'e_Toc',
        't2_75seqs':  'd_Toc',
        't3_27seqs':  'a_Toc',
        't4_92seqs':  'c_Toc',
        't5_31seqs':  'b2_Toc',
        't6_324seqs': 'b1_Toc',
    },
    'teu': {
        't1_45seqs':  'e_Teu',
        't2_75seqs':  'd_Teu',
        't3_27seqs':  'a_Teu',
        't4_92seqs':  'c_Teu',
        't5_31seqs':  'b2_Teu',
        't6_324seqs': 'b1_Teu',
    },
    'gpy': {
        'd1_16seqs':  'a_Gpy',
        'd2_8seqs':   'd_Gpy',
        'd3_38seqs':  'c_Gpy',
        'd4_266seqs': 'e_Gpy',
        'd5_268seqs': 'b_Gpy',
    },
    'dmo': {
        'd1_16seqs':  'a_Dmo',
        'd2_8seqs':   'd_Dmo',
        'd3_38seqs':  'c_Dmo',
        'd4_266seqs': 'e_Dmo',
        'd5_268seqs': 'b_Dmo',
    },
}

# ---------------------------------------------------------------------------
# 2. Color scheme  (fill = Table-1 exact hex; line = visible darkened version)
# ---------------------------------------------------------------------------
# Determine fill color from new name
TABLE1_FILL = {
    'a': '#FFC000',          # orange
    'b': '#FFFF00',          # yellow (b, b1, b2, b3)
}
# Species-specific for c–g
SP_FILL = {
    'toc': '#C6D9F1',   # light blue
    'teu': '#C6D9F1',
    'gpy': '#F2DBDB',   # light pink
    'dmo': '#F2DBDB',
    'saq': '#D6E3BC',   # light sage-green
    'ccr': '#E5DFEC',   # light lavender
}

# Dark line colors per hue group (for visibility on white background)
LINE_COLOR = {
    '#FFC000': '#B87000',   # darker amber
    '#FFFF00': '#909000',   # dark olive
    '#C6D9F1': '#1A5C96',   # steel blue
    '#F2DBDB': '#A93226',   # dark red
    '#D6E3BC': '#3A7D35',   # forest green
    '#E5DFEC': '#6A1FA0',   # violet
}

# Within-group line-color variants (for up to 5 subfamilies sharing one hue)
# Each group gets a graduated list from lighter to darker of the same hue
LINE_VARIANTS = {
    '#FFC000': ['#E6A000', '#CC8800', '#B87000', '#A06000', '#885000'],
    '#FFFF00': ['#B8B800', '#A0A000', '#909000', '#787800', '#606000'],
    '#C6D9F1': ['#5B9CD6', '#2E7CBD', '#1A5C96', '#0D3E6A', '#082647'],
    '#F2DBDB': ['#E07070', '#C05050', '#A93226', '#8B1A1A', '#6B0808'],
    '#D6E3BC': ['#8DC85A', '#5EAA3A', '#3A7D35', '#256020', '#16430F'],
    '#E5DFEC': ['#B07FD4', '#8B50BA', '#6A1FA0', '#4E0D80', '#380060'],
}


def get_fill_color(new_name, sp):
    letter = new_name[0].lower()   # 'a', 'b', 'c', 'd', ...
    if letter == 'a':
        return '#FFC000'
    if letter == 'b':
        return '#FFFF00'
    return SP_FILL.get(sp, '#888888')


def get_line_color(fill, variant_idx=0):
    variants = LINE_VARIANTS.get(fill, ['#444444'] * 5)
    return variants[min(variant_idx, len(variants) - 1)]


# ---------------------------------------------------------------------------
# 3. Parse Plotly.newPlot('plot_div_kde', ...) from HTML
# ---------------------------------------------------------------------------
def extract_matching(text, start, open_char, close_char):
    """Return substring from start that ends at matching close_char."""
    depth = 0
    i = start
    while i < len(text):
        if text[i] == open_char:
            depth += 1
        elif text[i] == close_char:
            depth -= 1
            if depth == 0:
                return text[start:i + 1]
        i += 1
    return None


def extract_kde_traces(html_text):
    """Extract traces list and layout dict from Plotly.newPlot('plot_div_kde',...) call."""
    marker = "Plotly.newPlot('plot_div_kde',"
    pos = html_text.find(marker)
    if pos == -1:
        return None, None
    pos += len(marker)
    # skip whitespace
    while html_text[pos] in ' \t\n\r':
        pos += 1
    # extract trace array
    if html_text[pos] != '[':
        return None, None
    traces_str = extract_matching(html_text, pos, '[', ']')
    if not traces_str:
        return None, None
    pos2 = pos + len(traces_str)
    # skip comma + whitespace to layout
    while html_text[pos2] in ', \t\n\r':
        pos2 += 1
    if html_text[pos2] != '{':
        return None, None
    layout_str = extract_matching(html_text, pos2, '{', '}')

    try:
        traces = json.loads(traces_str)
        layout = json.loads(layout_str)
    except json.JSONDecodeError as e:
        print(f"  JSON parse error: {e}")
        return None, None
    return traces, layout


# ---------------------------------------------------------------------------
# 4. Transform traces
# ---------------------------------------------------------------------------
def transform_traces(traces, sp):
    """Rename old codes, apply Table-1 colors, add fill."""
    rename_map = RENAME[sp]

    # group fills to count variants
    fill_counter = {}  # fill_color -> list of trace indices

    transformed = []
    for t in traces:
        if t.get('type') != 'scatter' or t.get('mode') != 'lines':
            continue
        old_name = t.get('name', '')
        new_name = rename_map.get(old_name, old_name)
        fill = get_fill_color(new_name, sp)
        fill_key = fill
        if fill_key not in fill_counter:
            fill_counter[fill_key] = []
        fill_counter[fill_key].append(len(transformed))

        nt = {
            'type': 'scatter',
            'mode': 'lines',
            'name': new_name,
            'x': t['x'],
            'y': t['y'],
            'fill': 'tozeroy',
            'fillcolor': _hex_to_rgba(fill, 0.25),
            'line': {'color': '#000000', 'width': 2},  # placeholder, updated below
            'hovertemplate': f'<b>{new_name}</b><br>divergence ~%{{x:.0f}}%<br>copies %{{y:,d}}<extra></extra>',
        }
        transformed.append((nt, fill))

    # Assign per-variant line colors within each fill group
    # First pass: count per fill
    fill_seencount = {}
    result = []
    for nt, fill in transformed:
        idx = fill_seencount.get(fill, 0)
        fill_seencount[fill] = idx + 1
        line_color = get_line_color(fill, idx)
        nt['line']['color'] = line_color
        # Legend group: cluster same-fill traces together in legend
        nt['legendgroup'] = fill
        result.append(nt)

    return result


def _hex_to_rgba(hex_color, alpha):
    h = hex_color.lstrip('#')
    r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    return f'rgba({r},{g},{b},{alpha})'


# ---------------------------------------------------------------------------
# 5. Species metadata
# ---------------------------------------------------------------------------
SPECIES_ORDER = ['toc', 'teu', 'gpy', 'dmo', 'saq', 'ccr']
PANEL_LABELS  = list('ABCDEF')

SPECIES_FULL = {
    'toc': 'Talpa occidentalis',
    'teu': 'Talpa europaea',
    'gpy': 'Galemys pyrenaicus',
    'dmo': 'Desmana moschata',
    'saq': 'Scalopus aquaticus',
    'ccr': 'Condylura cristata',
}

REPORT_PATHS = {sp: os.path.join(sp, 'report.html') for sp in SPECIES_ORDER}

# ---------------------------------------------------------------------------
# 6. Build legend HTML
# ---------------------------------------------------------------------------
LEGEND_GROUPS = [
    ('#FFC000', '#B87000', 'a subfamilies (oldest retained lineage, all species)',
     ['a_Toc', 'a_Teu', 'a_Gpy', 'a_Dmo', 'a_Saq']),
    ('#FFFF00', '#909000', 'b subfamilies (all species)',
     ['b1_Toc', 'b2_Toc', 'b1_Teu', 'b2_Teu', 'b_Gpy', 'b_Dmo',
      'b1_Saq', 'b2_Saq', 'b1_Ccr', 'b2_Ccr']),
    ('#C6D9F1', '#1A5C96', 'c–e subfamilies: Talpa spp. (toc, teu)',
     ['c_Toc', 'd_Toc', 'e_Toc', 'c_Teu', 'd_Teu', 'e_Teu']),
    ('#F2DBDB', '#A93226', 'c–e subfamilies: Galemys / Desmana (gpy, dmo)',
     ['c_Gpy', 'd_Gpy', 'e_Gpy', 'c_Dmo', 'd_Dmo', 'e_Dmo']),
    ('#D6E3BC', '#3A7D35', 'c–g subfamilies: Scalopus (saq)',
     ['c_Saq (s1)', 'c_Saq (s9)', 'd_Saq', 'e_Saq', 'f_Saq', 'g_Saq']),
    ('#E5DFEC', '#6A1FA0', 'c–g subfamilies: Condylura (ccr)',
     ['c_Ccr', 'd_Ccr', 'e_Ccr', 'f_Ccr', 'g_Ccr']),
]


def build_legend_html():
    parts = ['<div class="legend-container">',
             '<h3>Legend — Subfamily colour groups (Table 1)</h3>']
    for fill, line, label, members in LEGEND_GROUPS:
        parts.append(f'''
  <div class="legend-row">
    <span class="swatch" style="background:{fill};border:2px solid {line};"></span>
    <span class="legend-label"><strong>{htmllib.escape(label)}:</strong>
      {htmllib.escape(", ".join(members))}</span>
  </div>''')
    parts.append('</div>')
    return '\n'.join(parts)


# ---------------------------------------------------------------------------
# 7. Assemble HTML
# ---------------------------------------------------------------------------
def build_html(panels_data):
    """panels_data: list of (sp, panel_label, traces_json, layout) in order."""

    panels_html = []
    for sp, label, traces, layout in panels_data:
        full_name = SPECIES_FULL[sp]
        traces_json = json.dumps(traces)
        layout_json = json.dumps({
            **layout,
            'title': None,
            'xaxis': {'title': 'Divergence from consensus (%)', 'range': [0, 45]},
            'yaxis': {'title': 'Copies'},
            'legend': {'orientation': 'v', 'x': 1.02, 'y': 1,
                       'xanchor': 'left', 'yanchor': 'top',
                       'tracegroupgap': 2},
            'margin': {'t': 20, 'b': 60, 'l': 70, 'r': 180},
            'height': 320,
        })
        panels_html.append(f'''
  <div class="panel" id="panel-{sp}" data-sp="{sp}" data-label="{label}">
    <div class="panel-header">
      <span class="panel-letter">{label}</span>
      <span class="panel-title"><em>{htmllib.escape(full_name)}</em></span>
      <div class="panel-btns">
        <button class="move-btn" onclick="movePanel('{sp}', -1)" title="Move up">&#8679;</button>
        <button class="move-btn" onclick="movePanel('{sp}', +1)" title="Move down">&#8681;</button>
        <button class="dl-btn" onclick="dlPNG('{sp}')" title="Download PNG">&#128247; PNG</button>
        <button class="dl-btn" onclick="dlSVG('{sp}')" title="Download SVG">&#128196; SVG</button>
      </div>
    </div>
    <div class="plot-wrap" id="plot-{sp}"></div>
  </div>''')

    init_js = []
    for sp, label, traces, layout in panels_data:
        traces_json = json.dumps(traces)
        layout_json = json.dumps({
            'title': None,
            'xaxis': {'title': 'Divergence from consensus (%)', 'range': [0, 45]},
            'yaxis': {'title': 'Copies'},
            'legend': {'orientation': 'v', 'x': 1.02, 'y': 1,
                       'xanchor': 'left', 'yanchor': 'top',
                       'tracegroupgap': 2},
            'margin': {'t': 20, 'b': 60, 'l': 70, 'r': 180},
            'height': 320,
        })
        init_js.append(
            f"  Plotly.newPlot('plot-{sp}', {traces_json}, {layout_json}, "
            f"{{responsive:true, displaylogo:false, modeBarButtonsToRemove:['toImage']}});"
        )

    legend_html = build_legend_html()

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Tal SINE — divergence from consensus (6 species)</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
  body {{ font-family: Arial, sans-serif; max-width: 1100px; margin: 0 auto; padding: 20px;
         background: #fff; color: #222; }}
  h1 {{ font-size: 1.3em; margin-bottom: 4px; }}
  .subtitle {{ color: #555; font-size: 0.9em; margin-bottom: 20px; }}
  /* panels */
  #figure-container {{ display:flex; flex-direction:column; gap:8px; }}
  .panel {{ border: 1px solid #ddd; border-radius: 4px; padding: 6px 10px 2px 10px;
            background: #fafafa; }}
  .panel-header {{ display:flex; align-items:center; gap:10px; margin-bottom:4px; }}
  .panel-letter {{ font-size: 1.4em; font-weight: bold; min-width: 22px; color:#333; }}
  .panel-title {{ font-size: 0.95em; flex:1; }}
  .panel-btns {{ display:flex; gap:4px; }}
  .move-btn {{ background:#e8e8e8; border:1px solid #bbb; border-radius:3px;
               cursor:pointer; padding:1px 6px; font-size:1.1em; line-height:1.2; }}
  .move-btn:hover {{ background:#d0d0d0; }}
  .dl-btn {{ background:#d4e8ff; border:1px solid #88b; border-radius:3px;
             cursor:pointer; padding:2px 7px; font-size:0.8em; }}
  .dl-btn:hover {{ background:#b8d0f0; }}
  .plot-wrap {{ width:100%; }}
  /* legend */
  .legend-container {{ margin-top:28px; padding:14px 18px; border:1px solid #ddd;
                       border-radius:4px; background:#f7f7f7; }}
  .legend-container h3 {{ margin:0 0 10px 0; font-size:1em; }}
  .legend-row {{ display:flex; align-items:center; margin:5px 0; }}
  .swatch {{ display:inline-block; width:28px; height:16px; border-radius:3px;
             margin-right:10px; flex-shrink:0; }}
  .legend-label {{ font-size:0.85em; }}
  /* caption + methods */
  .caption-block, .methods-block {{
    margin-top:24px; padding:14px 18px; border:1px solid #ddd;
    border-radius:4px; background:#f7f7f7;
  }}
  .caption-block h3, .methods-block h3 {{ margin:0 0 8px 0; font-size:1em; }}
  .caption-block p, .methods-block p {{ font-size:0.88em; line-height:1.6; margin:0; }}
  /* reorder highlight */
  .panel.dragging {{ opacity:0.5; }}
  /* download-all button */
  #dl-all-btn {{ margin-top:10px; padding:6px 16px; background:#4A8C3F; color:#fff;
                 border:none; border-radius:4px; cursor:pointer; font-size:0.9em; }}
  #dl-all-btn:hover {{ background:#3A7D35; }}
</style>
</head>
<body>

<h1>Tal SINEs — divergence from subfamily consensus across six Talpidae species</h1>
<p class="subtitle">
  Panels A–F can be reordered using the &#8679;&#8681; buttons.
  Click &#128247;&nbsp;PNG or &#128196;&nbsp;SVG to export individual panels.
  <button id="dl-all-btn" onclick="dlAll()">&#128247; Download all panels as PNG</button>
</p>

<div id="figure-container">
{''.join(panels_html)}
</div>

{legend_html}

<div class="caption-block">
  <h3>Figure caption</h3>
  <p>
    <strong>Figure N.</strong> Divergence of Tal SINE copies from their subfamily consensus sequences
    in six Talpidae species. Each panel shows copy-count histograms (0.5 percentage-point bins)
    per subfamily. Divergence was computed as 100&nbsp;&minus;&nbsp;(<em>copy bitscore /
    self-bitscore</em>&nbsp;&times;&nbsp;100%), capped at 0 (Ssearch36-based alignment,
    SINEderella pipeline step 2). Colour groups follow the cross-species subfamily nomenclature
    in Table&nbsp;1: amber&nbsp;=&nbsp;a (oldest), yellow&nbsp;=&nbsp;b, blue&nbsp;=&nbsp;c–e
    for <em>Talpa</em> spp., pink&nbsp;=&nbsp;c–e for <em>Galemys</em>/<em>Desmana</em>,
    green&nbsp;=&nbsp;c–g for <em>Scalopus</em>, lavender&nbsp;=&nbsp;c–g for
    <em>Condylura</em>. Up to 10,000 copies per subfamily were used.
    (A)&nbsp;<em>Talpa occidentalis</em>; (B)&nbsp;<em>Talpa europaea</em>;
    (C)&nbsp;<em>Galemys pyrenaicus</em>; (D)&nbsp;<em>Desmana moschata</em>;
    (E)&nbsp;<em>Scalopus aquaticus</em>; (F)&nbsp;<em>Condylura cristata</em>.
  </p>
</div>

<div class="methods-block">
  <h3>Methods — divergence metric</h3>
  <p>
    Divergence from consensus for each assigned SINE copy was computed as
    100&nbsp;&minus;&nbsp;(<em>bitscore<sub>copy→consensus</sub></em> /
    <em>bitscore<sub>consensus→self</sub></em>&nbsp;&times;&nbsp;100%).
    Values were clamped to 0 where the copy bitscore exceeded the self-bitscore
    (full-length high-quality copies). Assignments were performed by the
    SINEderella pipeline (Toki, unpubl.) using Ssearch36 (Smith–Waterman local
    alignment) against the per-species SINE bank of consensus sequences derived
    from step-2 clustering. All assigned copies per subfamily were binned into
    0.5 percentage-point divergence bins; bar height = copy count (raw, not
    normalised). Subfamilies were unified into cross-species homology groups
    a–g using pairwise consensus alignment scores and phylogenetic context
    (Table&nbsp;1). Colours in the figure reflect these homology groups as
    listed in the legend above; line shading within the same hue distinguishes
    individual subfamilies sharing the same group colour.
  </p>
</div>

<script>
// ---- Plotly initialisation ----
{chr(10).join(init_js)}

// ---- Reorder panels ----
function movePanel(sp, dir) {{
  const container = document.getElementById('figure-container');
  const panels = Array.from(container.children);
  const el = document.getElementById('panel-' + sp);
  const idx = panels.indexOf(el);
  const newIdx = idx + dir;
  if (newIdx < 0 || newIdx >= panels.length) return;
  if (dir === -1) container.insertBefore(el, panels[newIdx]);
  else container.insertBefore(panels[newIdx], el);
  relabelPanels();
}}

function relabelPanels() {{
  const labels = 'ABCDEF';
  const panels = Array.from(document.getElementById('figure-container').children);
  panels.forEach((p, i) => {{
    const letter = labels[i] || String(i + 1);
    p.querySelector('.panel-letter').textContent = letter;
    p.dataset.label = letter;
  }});
}}

// ---- PNG / SVG download ----
function dlPNG(sp) {{
  const sp2full = {json.dumps({sp: SPECIES_FULL[sp].replace(' ','_') for sp in SPECIES_ORDER})};
  Plotly.downloadImage('plot-' + sp, {{
    format: 'png', width: 1400, height: 500, scale: 2,
    filename: sp2full[sp] + '_divergence'
  }});
}}

function dlSVG(sp) {{
  const sp2full = {json.dumps({sp: SPECIES_FULL[sp].replace(' ','_') for sp in SPECIES_ORDER})};
  Plotly.downloadImage('plot-' + sp, {{
    format: 'svg', width: 1400, height: 500,
    filename: sp2full[sp] + '_divergence'
  }});
}}

function dlAll() {{
  const sps = {json.dumps(SPECIES_ORDER)};
  sps.forEach((sp, i) => {{
    setTimeout(() => dlPNG(sp), i * 800);
  }});
}}
</script>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# 8. Main
# ---------------------------------------------------------------------------
def main():
    os.makedirs('figures', exist_ok=True)

    panels_data = []
    for sp in SPECIES_ORDER:
        path = REPORT_PATHS[sp]
        print(f'Reading {path} ...', end=' ')
        with open(path, encoding='utf-8') as fh:
            html_text = fh.read()

        traces, layout = extract_kde_traces(html_text)
        if traces is None:
            print(f'ERROR: could not extract plot_div_kde from {path}')
            continue
        print(f'{len(traces)} traces found')

        new_traces = transform_traces(traces, sp)
        label = PANEL_LABELS[SPECIES_ORDER.index(sp)]
        panels_data.append((sp, label, new_traces, layout or {}))

        # Save JSON for kaleido export
        json_path = f'figures/plot_data_{sp}.json'
        with open(json_path, 'w', encoding='utf-8') as fh:
            json.dump({'traces': new_traces, 'sp': sp,
                       'full_name': SPECIES_FULL[sp]}, fh)
        print(f'  -> saved {json_path}')

    out_path = os.path.join('figures', 'conservation_figure.html')
    html_content = build_html(panels_data)
    with open(out_path, 'w', encoding='utf-8') as fh:
        fh.write(html_content)
    print(f'\nWrote {out_path}  ({len(html_content)//1024} KB)')


if __name__ == '__main__':
    main()
