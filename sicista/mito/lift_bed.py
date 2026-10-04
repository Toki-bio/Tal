#!/usr/bin/env python3
"""Lift a BED file from one sequence to another through an alignment that contains both.

    python3 lift_bed.py in.bed aln.fa SRC_NAME DST_NAME [DST_CHROM] > out.bed

SRC_NAME / DST_NAME are matched against the first word of the FASTA headers in
aln.fa. Each BED interval (0-based, half-open on SRC) is mapped to the alignment
columns of its first and last base and then to the DST bases at or inside those
columns; an interval that lands on nothing in DST is dropped (reported on
stderr). Column 1 of the output is DST_CHROM (default DST_NAME); the other
columns are copied, with thickStart/thickEnd and blockSizes reset to the lifted
interval.
"""
import sys


def read_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path, encoding='utf-8', errors='replace'):
        line = line.rstrip('\n')
        if line.startswith('>'):
            if name is not None:
                seqs[name] = ''.join(buf)
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        seqs[name] = ''.join(buf)
    return seqs


def main():
    bed, aln, src, dst = sys.argv[1:5]
    dst_chrom = sys.argv[5] if len(sys.argv) > 5 else dst
    seqs = read_fasta(aln)
    s, d = seqs[src], seqs[dst]
    src_cols = [i for i, c in enumerate(s) if c not in '-.']           # src base -> column
    col_to_dst = {}                                                     # column -> dst base
    k = 0
    for i, c in enumerate(d):
        if c not in '-.':
            col_to_dst[i] = k
            k += 1
    dst_len = k

    def first_dst_at_or_after(col):
        while col < len(d) and col not in col_to_dst:
            col += 1
        return col_to_dst.get(col)

    def last_dst_at_or_before(col):
        while col >= 0 and col not in col_to_dst:
            col -= 1
        return col_to_dst.get(col)

    dropped = 0
    for line in open(bed, encoding='utf-8'):
        if not line.strip() or line[0] == '#' or line.startswith(('browser', 'track')):
            if line.startswith('track'):
                sys.stdout.write(line.replace('"', '"', 1).rstrip('\n') + ' lifted_to="%s"\n' % dst_chrom)
            continue
        f = line.rstrip('\n').split('\t')
        a, b = int(f[1]), int(f[2])
        if a >= len(src_cols):
            dropped += 1
            continue
        b = min(b, len(src_cols))
        da = first_dst_at_or_after(src_cols[a])
        db = last_dst_at_or_before(src_cols[b - 1])
        if da is None or db is None or db < da:
            dropped += 1
            sys.stderr.write('dropped %s (%d-%d): no %s bases in those columns\n' % (f[3] if len(f) > 3 else '?', a, b, dst))
            continue
        f[0], f[1], f[2] = dst_chrom, str(da), str(db + 1)
        if len(f) > 7:
            f[6], f[7] = str(da), str(db + 1)
        if len(f) > 11:
            f[9], f[10], f[11] = '1', str(db + 1 - da), '0'
        print('\t'.join(f))
    sys.stderr.write('%s: %d bases; %d intervals dropped\n' % (dst, dst_len, dropped))


if __name__ == '__main__':
    main()
