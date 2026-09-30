#!/usr/bin/env python3
"""Download every file in data_prep/download_manifest.tsv into its folder under data/raw/ (run from the repository root).
4 files in parallel, curl with retries and resume (-C -). Safe to re-run: finished files resume instantly.
Progress lines go to stdout and _logs/download_status.tsv."""
import csv, os, subprocess, time
from concurrent.futures import ThreadPoolExecutor

MANIFEST = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'download_manifest.tsv')
HERE = os.path.abspath('data/raw')   # raw-data root; manifest folders are relative to it
os.makedirs(HERE, exist_ok=True)
os.chdir(HERE)
os.makedirs('_logs', exist_ok=True)
rows = list(csv.DictReader(open(MANIFEST), delimiter='\t'))
# skip files already logged OK in a previous run; OpenGWAS rows need a token and are fetched separately
done = set()
if os.path.exists('_logs/download_status.tsv'):
    done = {l.split('\t')[2] for l in open('_logs/download_status.tsv') if l.split('\t')[1] == 'OK'}
rows = [r for r in rows if os.path.join(r['folder'], r['filename']) not in done and 'opengwas' not in r['url']]
# optional priority order: PRIORITY="GCST006867,GCST006979" puts matching ids first, in that order
pr = [x for x in os.environ.get('PRIORITY', '').split(',') if x]
rows.sort(key=lambda r: pr.index(r['id']) if r['id'] in pr else len(pr))
print('to download:', len(rows), [r['id'] for r in rows][:10], flush=True)


def get(r):
    os.makedirs(r['folder'], exist_ok=True)
    dest = os.path.join(r['folder'], r['filename'])
    log = os.path.join('_logs', r['filename'].replace('/', '_') + '.log')
    for attempt in range(1, 4):
        t0 = time.time()
        p = subprocess.run(['curl', '-fL', '--retry', '5', '--retry-delay', '15', '--connect-timeout', '60',
                            '-C', '-', '-sS', '-o', dest, r['url']], capture_output=True, text=True)
        open(log, 'a').write(f'attempt {attempt} rc={p.returncode} {p.stderr.strip()}\n')
        # rc 33 = server refused range on an already complete file; treat as done if file exists
        if p.returncode == 0 or (p.returncode == 33 and os.path.getsize(dest) > 0):
            mb = os.path.getsize(dest) / 1e6
            line = f"OK\t{dest}\t{mb:.1f} MB\t{time.time() - t0:.0f}s"
            break
        time.sleep(30)
    else:
        line = f"FAIL\t{dest}\trc={p.returncode}\t{p.stderr.strip()[:200]}"
    print(line, flush=True)
    with open('_logs/download_status.tsv', 'a') as f:
        f.write(time.strftime('%Y-%m-%d %H:%M:%S') + '\t' + line + '\n')
    return line


with ThreadPoolExecutor(int(os.environ.get('NWORKERS', '4'))) as ex:
    res = list(ex.map(get, rows))
print('ALL_DONE', sum(r.startswith('OK') for r in res), 'OK /', len(res), flush=True)
