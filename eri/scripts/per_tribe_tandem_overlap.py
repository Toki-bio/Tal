#!/usr/bin/env python3
import re, subprocess, glob, os

RUN = '/data/W/toki/Genomes/Mammalia/Eulipotyphla/Erniacidae/run_20260820_221537'
TRIBES_DIR = f'{RUN}/eri_tribes_out'
TANDEM_BED = f'{RUN}/trf_out/tandem_repeats.merged.underscored.bed'

files = sorted(glob.glob(f'{TRIBES_DIR}/eri_tribe*.aln.fa'))
for fn in files:
    tribe = os.path.basename(fn).split('_')[1]
    bed_lines = []
    name = None
    seq = ''
    with open(fn) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                name = line[1:]
            # only need headers
    with open(fn) as f:
        for line in f:
            if line.startswith('>'):
                hdr = line[1:].strip()
                if 'CONSENSUS' in hdr:
                    continue
                m = re.match(r'^([^:]+):(\d+)-(\d+)\(([+-])\)', hdr)
                if not m:
                    continue
                ctg, s, e, strand = m.groups()
                bed_lines.append(f'{ctg}\t{s}\t{e}\t{hdr}\t0\t{strand}')
    bed_path = f'/tmp/{tribe}.bed'
    with open(bed_path, 'w') as f:
        f.write('\n'.join(bed_lines) + '\n')
    total = len(bed_lines)
    sorted_bed = f'/tmp/{tribe}.sorted.bed'
    subprocess.run(f'bedtools sort -i {bed_path} > {sorted_bed}', shell=True, check=True)
    result = subprocess.run(f'bedtools intersect -u -a {sorted_bed} -b {TANDEM_BED} | wc -l',
                             shell=True, capture_output=True, text=True, check=True)
    tand = int(result.stdout.strip())
    pct = 100*tand/total if total else 0
    print(f'{tribe}: total={total} tandem={tand} pct={pct:.1f}%')
