"""
Export each species divergence panel as PNG (+ TIFF if PIL available).

Requirements: pip install plotly kaleido Pillow

Run from c:\work\Tal:
    .\.venv\Scripts\python.exe _export_divergence_pngs.py
"""

import json, os, sys

try:
    import plotly.graph_objects as go
    import plotly.io as pio
except ImportError:
    sys.exit('plotly not found — run: pip install plotly kaleido')

try:
    import kaleido  # noqa: F401 — just check it's installed
except ImportError:
    sys.exit('kaleido not found — run: pip install kaleido')

os.makedirs('figures', exist_ok=True)

LAYOUT_BASE = dict(
    xaxis=dict(title='Divergence from consensus (%)', range=[0, 45]),
    yaxis=dict(title='Copies'),
    legend=dict(orientation='v', x=1.02, y=1,
                xanchor='left', yanchor='top', tracegroupgap=2),
    margin=dict(t=40, b=70, l=80, r=220),
    paper_bgcolor='white',
    plot_bgcolor='white',
)

SPECIES_FULL = {
    'toc': 'Talpa occidentalis',
    'teu': 'Talpa europaea',
    'gpy': 'Galemys pyrenaicus',
    'dmo': 'Desmana moschata',
    'saq': 'Scalopus aquaticus',
    'ccr': 'Condylura cristata',
}
LABELS = list('ABCDEF')
SPECIES_ORDER = ['toc', 'teu', 'gpy', 'dmo', 'saq', 'ccr']

for i, sp in enumerate(SPECIES_ORDER):
    json_path = f'figures/plot_data_{sp}.json'
    if not os.path.exists(json_path):
        print(f'SKIP {sp}: {json_path} not found (run _make_figure_page.py first)')
        continue

    with open(json_path) as fh:
        data = json.load(fh)

    traces = []
    for t in data['traces']:
        traces.append(go.Scatter(
            x=t['x'], y=t['y'],
            name=t['name'],
            mode='lines',
            fill=t.get('fill', 'none'),
            fillcolor=t.get('fillcolor', 'rgba(0,0,0,0.1)'),
            line=dict(color=t['line']['color'], width=t['line'].get('width', 2)),
            legendgroup=t.get('legendgroup', ''),
            hovertemplate=t.get('hovertemplate', ''),
        ))

    full_name = SPECIES_FULL[sp]
    label = LABELS[i]
    layout = go.Layout(
        title=dict(text=f'({label}) <i>{full_name}</i>',
                   font=dict(size=14), x=0.02, xanchor='left'),
        **LAYOUT_BASE,
    )

    fig = go.Figure(data=traces, layout=layout)

    safe_name = full_name.replace(' ', '_')

    # PNG (width=1400, height=500, scale=2 → 2800×1000 px, ~300 DPI for 9.3×3.3 in)
    png_path = f'figures/{safe_name}_divergence.png'
    print(f'Writing {png_path} ...', end=' ')
    try:
        pio.write_image(fig, png_path, format='png', width=1400, height=500, scale=2)
        print('OK')
    except Exception as e:
        print(f'ERROR: {e}')

    # TIFF via PIL (convert from PNG bytes)
    try:
        from PIL import Image
        import io
        png_bytes = pio.to_image(fig, format='png', width=1400, height=500, scale=2)
        img = Image.open(io.BytesIO(png_bytes))
        tiff_path = f'figures/{safe_name}_divergence.tiff'
        img.save(tiff_path, compression='tiff_lzw', dpi=(300, 300))
        print(f'  TIFF: {tiff_path}  OK')
    except ImportError:
        print(f'  (PIL not available — skipping TIFF for {sp})')
    except Exception as e:
        print(f'  TIFF ERROR: {e}')

print('\nDone.')
