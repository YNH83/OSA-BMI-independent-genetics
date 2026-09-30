"""Compare standardized files before (data/hm3_std) and after (data/hm3_std_fixP) the P-underflow fix (deviation record 03).
Run from the repository root."""
import json, glob, os
print(f"{'trait':26s} {'kept_old':>9s} {'kept_new':>9s} {'recov':>6s} {'gws_old':>7s} {'gws_new':>7s}")
for f in sorted(glob.glob('data/hm3_std/*.meta.json')):
    t = os.path.basename(f)[:-10]; g = f'data/hm3_std_fixP/{t}.meta.json'
    if not os.path.exists(g): continue
    a = json.load(open(f)); b = json.load(open(g))
    flag = '  <== CHANGED' if a['n_hm3_kept'] != b['n_hm3_kept'] else ''
    print(f"{t:26s} {a['n_hm3_kept']:9d} {b['n_hm3_kept']:9d} {b.get('n_p_underflow_recovered',0):6d} {a['n_gws']:7d} {b['n_gws']:7d}{flag}")
