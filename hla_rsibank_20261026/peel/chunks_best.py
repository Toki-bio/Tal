"""Best bank-or-owner consensus for every chunk consensus of the hla 30k SubFam plate, in plate order.
Usage: chunks_best.py <chunks.fa> <bank_plus_owner.fa> <chunks_vs_bank.m8> <plate_order.tsv> > chunks_best.tsv
Columns: plate row, chunk, chunk length (N removed), best consensus, % identity, aligned span on the chunk and on the
consensus, consensus length, second-best consensus and its bitscore as a fraction of the best."""
import collections, sys


def rd(p):
    d, n = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"):
            n = l[1:].split()[0]; d[n] = ""
        else:
            d[n] += l
    return d


ch, bk = rd(sys.argv[1]), rd(sys.argv[2])
hits = collections.defaultdict(list)
for l in open(sys.argv[3]):
    f = l.split("\t")
    hits[f[0]].append((float(f[11]), f[1], float(f[2]), int(f[6]), int(f[7]), int(f[8]), int(f[9])))
print("row\tchunk\tlen\tbest\tpid\tchunk_span\tcons_span\tcons_len\tsecond\tsecond_bits_ratio")
for l in open(sys.argv[4]):
    r, n = l.rstrip("\n").split("\t")
    if n not in ch:
        continue
    best = {}
    for h in sorted(hits.get(n, []), reverse=True):
        best.setdefault(h[1], h)
    B = sorted(best.values(), reverse=True)
    if not B:
        print("%s\t%s\t%d\t-" % (r, n, len(ch[n]))); continue
    b, s = B[0], (B[1] if len(B) > 1 else None)
    print("%s\t%s\t%d\t%s\t%.1f\t%d-%d\t%d-%d\t%d\t%s\t%s" % (r, n, len(ch[n]), b[1], b[2], b[3], b[4], b[5], b[6], len(bk[b[1]]),
                                                         s[1] if s else "-", "%.2f" % (s[0] / b[0]) if s else "-"))
