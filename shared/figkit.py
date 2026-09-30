"""Shared figure and raw-data helpers.

Palette: dataviz reference categorical slots 1 to 4 (validated 2026-09-26, light mode: CVD worst adjacent
dE 9.1, normal-vision 22.9; contrast WARN for aqua and yellow -> every figure carries direct labels).
Every display item is saved as 600 dpi RGB TIFF + PDF, and gets a raw-data Excel with sheets
data / sources / columns / checks.
"""
import hashlib, os, time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

BLUE, ORANGE, AQUA, YELLOW = '#2a78d6', '#eb6834', '#1baf7a', '#eda100'
GREY, INK, INK2 = '#8a8984', '#0b0b0b', '#52514e'
# Component colours (design spec rule 1): raw target trait (OSA) = INK2; the subtracted trait's own GWAS
# (BMI) = ORANGE_LIGHT; the component through it (F_BMI) = ORANGE; the residual (F_NonBMI) = BLUE.
ORANGE_LIGHT = '#f2a07b'
RG_POS, RG_NEG = '#9e2a55', '#3f7f55'   # rg diverging scale (wine to green, white midpoint); never blue/red (rule 3)
SEQ = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b']
DIV_NEG, DIV_MID, DIV_POS = '#2a78d6', '#f0efec', '#e34948'

plt.rcParams.update({
    'font.family': 'Arial', 'font.size': 8, 'axes.titlesize': 9, 'axes.labelsize': 8,
    'axes.edgecolor': INK2, 'axes.linewidth': 0.6, 'xtick.color': INK2, 'ytick.color': INK2,
    'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 300, 'pdf.fonttype': 42,
    'legend.frameon': False,
})


def sha256(path, n=None):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    d = h.hexdigest()
    return d[:n] if n else d


def save_fig(fig, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    paths = []
    for ext in ('tiff', 'pdf'):
        p = os.path.join(outdir, f'{name}.{ext}')
        kw = dict(pil_kwargs={'compression': 'tiff_lzw'}) if ext == 'tiff' else {}
        fig.savefig(p, dpi=600, bbox_inches='tight', **kw)
        if ext == 'tiff':  # submission systems reject alpha channels: store RGB at 600 dpi (see export_tiff600.py)
            from PIL import Image
            Image.open(p).convert('RGB').save(p, compression='tiff_lzw', dpi=(600, 600))
        paths.append(p)
    plt.close(fig)
    return paths


def write_raw_xlsx(path, data, inputs, script, columns, checks=None, notes=None):
    """data: DataFrame or {sheet_name: DataFrame}; inputs: list of file paths; script: producing script path;
    columns: {column: definition}; checks: list of (description, value) a reviewer can recompute by hand."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    sheets = data if isinstance(data, dict) else {'data': data}
    src = [dict(role='input', path=p, sha256=sha256(p) if os.path.exists(p) else 'MISSING') for p in inputs]
    src.append(dict(role='script', path=script, sha256=sha256(script)))
    src.append(dict(role='generated', path=path, sha256=time.strftime('%Y-%m-%d %H:%M:%S')))
    with pd.ExcelWriter(path, engine='openpyxl') as xw:
        for name, df in sheets.items():
            df.to_excel(xw, sheet_name=name[:31], index=False)
            for row in xw.book[name[:31]].iter_rows():   # strings such as '=~' (lavaan operator) must stay text, not formulas
                for c in row:
                    if isinstance(c.value, str) and c.value.startswith('='): c.data_type = 's'
        pd.DataFrame(src).to_excel(xw, sheet_name='sources', index=False)
        pd.DataFrame([dict(column=k, definition=v) for k, v in columns.items()]).to_excel(xw, sheet_name='columns', index=False)
        pd.DataFrame(checks or [], columns=['check', 'value']).to_excel(xw, sheet_name='checks', index=False)
        if notes:
            pd.DataFrame({'notes': notes}).to_excel(xw, sheet_name='notes', index=False)
    return path
