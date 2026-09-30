#!/usr/bin/env python3
"""Segmented, resumable download for servers that throttle per connection (EBI FTP, Zenodo).

Keeps any existing partial file as a valid prefix, splits the remaining bytes into N ranges fetched in
parallel (each range resumable in its own .partK file), then appends the parts in order and checks the
final size against the server's Content-Length.

Usage:
  seg_download.py --manifest-ids GCST006867 GCST009760 ...   (files from data_prep/download_manifest.tsv, in this order)
  seg_download.py --url URL --dest PATH                       (single file, e.g. reference data)
Env: NSEG (default 8).
"""
import argparse, csv, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

MANIFEST = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'download_manifest.tsv')
HERE = os.path.abspath('data/raw')   # raw-data root (run from the repository root)
NSEG = int(os.environ.get('NSEG', '8'))


def total_size(url):
    h = subprocess.run(['curl', '-sIL', '-m', '120', url], capture_output=True, text=True).stdout
    v = [l.split(':', 1)[1].strip() for l in h.splitlines() if l.lower().startswith('content-length')]
    return int(v[-1]) if v else None


def fetch_range(url, part, start, end):
    need = end - start + 1
    for attempt in range(30):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if have == need:
            return True
        if have > need:
            os.remove(part); have = 0
        with open(part, 'ab') as f:
            subprocess.run(['curl', '-sL', '--connect-timeout', '60', '-m', '1800', '-r', f'{start + have}-{end}', url],
                           stdout=f, stderr=subprocess.DEVNULL)
        time.sleep(3)
    return os.path.getsize(part) == need


def download(url, dest):
    t0 = time.time()
    os.makedirs(os.path.dirname(dest) or '.', exist_ok=True)
    tot = total_size(url)
    if not tot:
        return f'FAIL\t{dest}\tno content-length'
    have = os.path.getsize(dest) if os.path.exists(dest) else 0
    if have > tot:
        return f'FAIL\t{dest}\tlocal {have} > remote {tot}'
    if have < tot:
        rem = tot - have
        step = -(-rem // NSEG)
        ranges = [(have + i * step, min(have + (i + 1) * step - 1, tot - 1)) for i in range(NSEG) if have + i * step < tot]
        parts = [f'{dest}.part{i}' for i in range(len(ranges))]
        with ThreadPoolExecutor(len(ranges)) as ex:
            ok = list(ex.map(lambda a: fetch_range(url, *a), [(parts[i], *ranges[i]) for i in range(len(ranges))]))
        if not all(ok):
            return f'FAIL\t{dest}\tsegments incomplete (kept for resume)'
        with open(dest, 'ab') as out:
            for p in parts:
                with open(p, 'rb') as f:
                    while True:
                        b = f.read(1 << 22)
                        if not b: break
                        out.write(b)
        for p in parts: os.remove(p)
    size = os.path.getsize(dest)
    status = 'OK' if size == tot else 'FAIL'
    return f'{status}\t{dest}\t{size / 1e6:.1f} MB\t{time.time() - t0:.0f}s\t(seg {NSEG}, remote {tot})'


def log(line):
    print(line, flush=True)
    os.makedirs(os.path.join(HERE, '_logs'), exist_ok=True)
    with open(os.path.join(HERE, '_logs', 'download_status.tsv'), 'a') as f:
        f.write(time.strftime('%Y-%m-%d %H:%M:%S') + '\t' + line + '\n')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest-ids', nargs='*')
    ap.add_argument('--url'); ap.add_argument('--dest')
    a = ap.parse_args()
    if a.url:
        log(download(a.url, a.dest))
        sys.exit()
    rows = list(csv.DictReader(open(MANIFEST), delimiter='\t'))
    for gid in a.manifest_ids:
        for r in rows:
            if r['id'] == gid:
                log(download(r['url'], os.path.join(HERE, r['folder'], r['filename'])))
