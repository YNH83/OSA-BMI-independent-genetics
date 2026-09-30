#!/usr/bin/env python3
"""Verify downloaded files: full gzip decompression test, first data line, row count for gz text,
SHA256, and md5 cross-check where the source publishes md5sum.txt. Writes data/raw/SHA256SUMS.txt and data/raw/_logs/verify_report.tsv.
Run from the repository root."""
import csv, gzip, hashlib, os, zipfile, tarfile

MANIFEST = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'download_manifest.tsv')
HERE = os.path.abspath('data/raw')   # raw-data root; manifest folders are relative to it
os.chdir(HERE)
os.makedirs('_logs', exist_ok=True)
rows = list(csv.DictReader(open(MANIFEST), delimiter='\t'))
rep, sums = [], []
for r in rows:
    p = os.path.join(r['folder'], r['filename'])
    if not os.path.exists(p):
        rep.append((p, 'MISSING', '', '', '')); continue
    size = os.path.getsize(p)
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            sha.update(b); md5.update(b)
    status, lines, first = 'OK', '', ''
    try:
        if p.endswith('.tar.gz'):
            with tarfile.open(p) as t:
                lines = ';'.join(f'{m.name}:{m.size}' for m in t.getmembers() if m.isfile())
        elif p.endswith('.gz') or p.endswith('.tbl.gz'):
            n = 0
            with gzip.open(p, 'rt', errors='replace') as g:
                for i, line in enumerate(g):
                    if i == 1: first = line[:120].rstrip()
                    n += 1
            lines = str(n)
        elif p.endswith('.zip'):
            with zipfile.ZipFile(p) as z:
                bad = z.testzip()
                if bad: status = 'CORRUPT ' + bad
                lines = ';'.join(f'{i.filename}:{i.file_size}' for i in z.infolist())
    except Exception as e:
        status = 'CORRUPT ' + str(e)[:80]
    rep.append((p, status, f'{size / 1e6:.1f} MB', lines, first))
    sums.append(f'{sha.hexdigest()}  {p}')
    r['_md5'] = md5.hexdigest()
    print(status, p, f'{size / 1e6:.1f} MB', lines[:60], flush=True)

# md5 cross-check against publisher md5sum.txt files
for r in rows:
    folder = r['folder']; ref = os.path.join(folder, 'md5sum.txt')
    if r['filename'] == 'md5sum.txt' or not os.path.exists(ref) or '_md5' not in r: continue
    pub = {l.split()[-1].lstrip('*').split('/')[-1]: l.split()[0] for l in open(ref) if l.strip()}
    if r['filename'] in pub:
        ok = pub[r['filename']] == r['_md5']
        print('MD5', 'MATCH' if ok else 'MISMATCH', r['filename'])
        rep.append((os.path.join(folder, r['filename']), 'MD5 ' + ('MATCH' if ok else 'MISMATCH'), '', '', ''))

open('SHA256SUMS.txt', 'w').write('\n'.join(sums) + '\n')
with open('_logs/verify_report.tsv', 'w') as f:
    f.write('file\tstatus\tsize\trows_or_members\tfirst_data_line\n')
    for x in rep: f.write('\t'.join(x) + '\n')
print('DONE', sum(1 for x in rep if x[1] == 'OK'), 'OK of', len(rows), 'files')
