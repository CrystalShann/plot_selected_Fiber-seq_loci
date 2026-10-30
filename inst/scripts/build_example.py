#!/usr/bin/env python3
"""Deterministic complete-record subset. Requires tabix/bgzip on PATH."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess

p = argparse.ArgumentParser()
p.add_argument('--sample-sheet', required=True)
p.add_argument('--tss-bed', required=True)
p.add_argument('--tf-root', required=True)
p.add_argument('--output', required=True)
p.add_argument('--tabix', default='tabix')
p.add_argument('--bgzip', default='bgzip')
a = p.parse_args()
out = Path(a.output)
out.mkdir(parents=True, exist_ok=True)
with open(a.sample_sheet) as f:
    samples = list(csv.DictReader(f, delimiter='\t'))[:2]
regions = []
with open(a.tss_bed) as f:
    for line in f:
        c = line.rstrip().split('\t')
        if len(c) >= 6 and len(c[3].split(';')) >= 3 and c[3].split(';')[2] in ('ALAS1', 'PRDX1'):
            regions.append(dict(region_id=c[3].split(';')[2], chr=c[0], start=int(c[1])-1000,
                end=int(c[1])+1000, strand=c[5], tss=int(c[1])+1,
                analysis_start=int(c[1])-999, analysis_end=int(c[1])+1000, width=2000))
regions.sort(key=lambda r: r['region_id'])
assert len(regions) == 2

def query(path, r):
    text = subprocess.check_output([a.tabix, str(path), f"{r['chr']}:{r['start']+1}-{r['end']}"], universal_newlines=True)
    return [x.split('\t') for x in text.splitlines() if x]

provenance = []
sheet = []
for s in samples:
    name = s['sample_name']; root = Path(s['fire_dir'])
    row = dict(sample_name=name, sample_label=s['sample_label'])
    files = {}
    for r in regions:
        chrom = r['chr']
        sources = dict(m6a=root/'extracted_results'/'m6a_by_chr'/f'{name}.ft_extracted_m6a.{chrom}.bed.gz',
            nucleosome=root/'extracted_results'/'nuc_by_chr'/f'{name}.ft_extracted_nuc.{chrom}.bed.gz',
            tf=Path(a.tf_root)/name/f'{name}.recalled_tf.bed.gz',
            fire=root/'additional-outputs-v0.1'/'fire-peaks'/f'{name}-v0.1-fire-elements.bed.gz')
        m6a = query(sources['m6a'], r)
        eligible = sorted({x[3] for x in m6a if int(x[1]) <= r['start'] and int(x[2]) >= r['end']})[:100]
        assert eligible, f'No spanning reads: {name} {r}'
        for kind, path in sources.items():
            records = m6a if kind == 'm6a' else query(path, r)
            kept = [x for x in records if x[3] in eligible]
            files.setdefault((kind, chrom), set()).update('\t'.join(x) for x in kept)
            provenance.append(dict(sample=name, region=r['region_id'], track=kind,
                source=str(path), selected_read_ids=eligible, record_count=len(kept)))
    for (kind, chrom), records in files.items():
        dest = out/'tracks'/name/f'{kind}.{chrom}.bed'
        dest.parent.mkdir(parents=True, exist_ok=True)
        ordered = sorted(records, key=lambda x: (x.split('\t')[0], int(x.split('\t')[1]), int(x.split('\t')[2]), x))
        dest.write_text('\n'.join(ordered)+'\n')
        subprocess.run([a.bgzip, '-f', str(dest)], check=True)
        subprocess.run([a.tabix, '-f', '-p', 'bed', str(dest)+'.gz'], check=True)
        row[kind] = f'tracks/{name}/{kind}.{{chr}}.bed.gz'
        row[kind+'_format'] = dict(m6a='bed12_fibertools', nucleosome='bed12_fibertools', tf='bed15_fiberhmm_tf', fire='bed11_fire')[kind]
    sheet.append(row)
for filename, rows in [('samples.tsv', sheet), ('regions.tsv', regions)]:
    with open(out/filename, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
(out/'provenance.json').write_text(json.dumps(dict(selection='First two sample-sheet rows; lexicographically first <=100 full-span read IDs per sample and locus; complete overlapping source records retained', records=provenance), indent=2)+'\n')
(out/'SHA256SUMS').write_text(''.join(f'{hashlib.sha256(f.read_bytes()).hexdigest()}  {f.relative_to(out)}\n' for f in sorted(out.rglob('*.gz*'))))
