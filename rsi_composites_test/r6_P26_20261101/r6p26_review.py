"""r6p26_review.py RUN FS9_STATUS OUTDIR - review plate for r6_210seqs vs r5_r6_P26 (P26 = r5-type head 2-137 + r6 138-376).

Groups (step 2 of RUN, vote_summary + assignment_full):
  A  firm r6        : assigned r6_210seqs, 10/10 votes
  B  split r6/P26   : votes for both r6 and P26 in the 10 cycles (not assigned: no unanimous vote)
  C  P26, runner r6 : assigned r5_r6_P26 10/10, runner-up r6_210seqs
  D  firm P26       : assigned r5_r6_P26 10/10, runner-up not r6 (reference: P26 copies the vote never doubted)
Only copies whose flanks fs9 calls unique (FS9_STATUS: copy TAB status) are shown; up to N per group, random (seed 1).

Per copy: the genome window locus + 300 bp upstream + 100 bp downstream (copy orientation) is aligned to the P26 consensus
(ssearch36, local). Row: P26 columns 1-376 hold the copy's bases at the P26 positions they align to (upper case; insertions in the copy left
out, counted in ins=); the unaligned rest of the window in lower case, packed against the alignment's ends (upstream ends
right before the first aligned column, downstream starts right after the last). Rows sorted by group, then first P26 column.
Head test: the window's sequence upstream of its r6 part (before the base aligned to P26 138) against P26 1-137 alone:
head=<aligned bp>@<identity>%/<bits>; for comparison every copy's upstream is also scored against 20 dinucleotide
shuffles of P26 1-137 (head_null = best bits of the shuffles).
Out: OUTDIR/r6_P26_review.aln.fa, OUTDIR/r6_P26_review.tsv (one row per shown copy), OUTDIR/summary.tsv"""
import collections, csv, gzip, os, random, re, subprocess, sys

run, fs9st, out = sys.argv[1:4]
N = int(os.environ.get("N_PER_GROUP", "40"))
UP, DOWN = 300, 100
S2 = os.path.join(run, "step2", "step2_output")
os.makedirs(out, exist_ok=True)
tmp = os.path.join(out, "tmp"); os.makedirs(tmp, exist_ok=True)
R6, P26 = "r6_210seqs", "r5_r6_P26"


def fasta(p):
    s, n = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"): n = l[1:].split()[0]; s[n] = []
        elif n: s[n].append(l)
    return {k: "".join(v).upper() for k, v in s.items()}


cons = fasta(os.path.join(run, "consensuses.clean.fa"))
votes = collections.defaultdict(dict)
for l in gzip.open(os.path.join(S2, "vote_summary.tsv.gz"), "rt"):
    f = l.rstrip("\n").split("\t")
    votes[f[0]][f[1]] = int(f[2])
af = {}
for r in csv.DictReader(open(os.path.join(S2, "assignment_full.tsv")), delimiter="\t"):
    runner = (re.search(r"runner=([^;]+)", r["RunnerInfo"] or "") or [None, "NA"])[1]
    af[r["Sequence"]] = (r["Subfamily"], r["Status"], int(r["Votes"]), runner)
status = dict(l.rstrip("\n").split("\t")[:2] for l in open(fs9st))

groups = collections.defaultdict(list)
for c, (sf, st, v, runner) in af.items():
    vc = votes.get(c, {})
    if st == "assigned" and sf == R6 and v == 10: g = "A"
    elif R6 in vc and P26 in vc: g = "B"
    elif st == "assigned" and sf == P26 and v == 10: g = "C" if runner == R6 else "D"
    else: continue
    groups[g].append(c)
rng = random.Random(1)
shown, counts = {}, {}
for g in "ABCD":
    cs = sorted(groups[g]); uniq = [c for c in cs if status.get(c) == "unique"]
    counts[g] = (len(cs), len(uniq), collections.Counter(status.get(c, "not_tested") for c in cs))
    shown[g] = rng.sample(uniq, min(N, len(uniq)))


def locus(n):
    m = re.match(r"^(.+):(\d+)-(\d+)\(([+-])\)$", n)
    return m.group(1), int(m.group(2)), int(m.group(3)), m.group(4)


# windows: locus + UP upstream + DOWN downstream in copy orientation; record where the locus sits inside the window
bed = os.path.join(tmp, "win.bed"); fai = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(os.path.join(run, "genome.clean.fa.fai"))}
inner = {}
with open(bed, "w") as fo:
    for g in "ABCD":
        for c in shown[g]:
            ctg, s, e, st = locus(c)          # names are 0-based start like BED (step 1 names); checked against getfasta below
            a, b = (s - UP, e + DOWN) if st == "+" else (s - DOWN, e + UP)
            a, b = max(0, a), min(fai[ctg], b)
            lo = (s - a) if st == "+" else (b - e)          # bases before the locus in copy orientation
            inner[c] = (lo, lo + (e - s))
            fo.write("%s\t%d\t%d\t%s\t0\t%s\n" % (ctg, a, b, c, st))
subprocess.run("bedtools getfasta -s -nameOnly -fi %s -bed %s > %s" % (os.path.join(run, "genome.clean.fa"), bed, os.path.join(tmp, "win.fa")), shell=True, check=True)
win = {re.sub(r"\([+-]\)$", "", k): v for k, v in fasta(os.path.join(tmp, "win.fa")).items()}


def markx10(q, qs, t, ts):
    open(os.path.join(tmp, "q.fa"), "w").write(">%s\n%s\n" % (q, qs)); open(os.path.join(tmp, "t.fa"), "w").write(">%s\n%s\n" % (t, ts))
    o = subprocess.run(["ssearch36", "-3", "-m", "10", "-b", "1", "-d", "1", os.path.join(tmp, "q.fa"), os.path.join(tmp, "t.fa")], capture_output=True, text=True, check=True).stdout
    if "\n>>" not in o: return None
    part = o.split("\n>>", 1)[1].split("\n>>>", 1)[0]
    bits = float(re.search(r"; sw_bits: (\S+)", part).group(1)); ident = float(re.search(r"; sw_ident: (\S+)", part).group(1))
    recs = []
    for blk in re.split(r"\n>(?!>)", part)[1:3]:
        blk = blk.split("\n; al_cons")[0]
        info = dict(re.findall(r"^; (\w+): (\S+)", blk, re.M))
        disp = "".join(x.strip() for x in blk.split("\n")[1:] if x.strip() and not x.startswith(";"))
        recs.append((int(info["al_start"]), int(info["al_stop"]), int(info["al_display_start"]), disp))
    return bits, ident, recs


def pairs(recs):
    """aligned (query pos, target pos) pairs inside the local alignment, 1-based; query insertions as (qpos, None)"""
    (qa, qb, qd, qdisp), (ta, tb, td, tdisp) = recs
    qk, tk, out = qd - 1, td - 1, []
    for j in range(max(len(qdisp), len(tdisp))):
        x = qdisp[j] if j < len(qdisp) else "-"; y = tdisp[j] if j < len(tdisp) else "-"
        if x != "-": qk += 1
        if y != "-": tk += 1
        if qa <= qk <= qb and ta <= tk <= tb and x != "-":
            out.append((qk, tk if y != "-" else None))
    return out


P = cons[P26]; L = len(P); HEAD = P[:137]
rnd = random.Random(7)
def dishuf(s):          # dinucleotide-preserving shuffle (Altschul-Erickson by Euler path is heavy; swap-based approx)
    s = list(s)
    for _ in range(len(s) * 10):
        i, j = sorted(rnd.sample(range(1, len(s) - 2), 2))
        if j - i > 1 and s[i - 1] == s[j - 1] and s[i + 1] == s[j + 1]: s[i], s[j] = s[j], s[i]
    return "".join(s)
shufs = [dishuf(HEAD) for _ in range(20)]

rows, info = [], []
for g in "ABCD":
    for c in shown[g]:
        w = win.get(c)
        if not w: continue
        lo, hi = inner[c]
        al = markx10(c, w, P26, P)
        if not al: continue
        bits, ident, recs = al
        pr = pairs(recs)
        qmap = {t: q for q, t in pr if t}
        ins = sum(1 for q, t in pr if t is None)
        qs_, qe_ = recs[0][0], recs[0][1]
        tpos = sorted(qmap)
        # head test: upstream of the base aligned to P26 138 (or of the alignment start if 138 is not aligned)
        cut = qmap.get(138) or (qs_ if tpos and tpos[0] >= 138 else None)
        up = w[:cut - 1] if cut else w[:lo]
        h = markx10("up", up, "head", HEAD) if len(up) >= 20 else None
        hb, hi_, hl = (h[0], h[1] * 100, h[2][1][1] - h[2][1][0] + 1) if h else (0.0, 0.0, 0)
        hnull = max((markx10("up", up, "shuf", sh) or (0.0,))[0] for sh in shufs) if len(up) >= 20 else 0.0
        # aligned to P26 = upper case at its P26 column; unaligned = lower case, packed against the alignment's ends
        # (upstream ends right before the first aligned column, downstream starts right after the last), as on the
        # flankscan plates: a copy without the r5 head shows its own flank in the head columns, in lower case
        cell = ["-"] * L
        for t, q in qmap.items(): cell[t - 1] = w[q - 1].upper()
        upstr, downstr = w[:qs_ - 1].lower(), w[qe_:].lower()
        vc = votes.get(c, {})
        vtxt = "r6:%d_P26:%d" % (vc.get(R6, 0), vc.get(P26, 0))
        name = "%s|%s|%s|P26span=%d-%d|head=%dbp@%.0f%%/%.0fb_null%.0fb|ins=%d" % (g, vtxt, c.replace("@U@", "_"), tpos[0], tpos[-1], hl, hi_, hb, hnull, ins)
        rows.append((g, tpos[0], name, upstr, "".join(cell), downstr, tpos[0], tpos[-1]))
        info.append((g, c, vtxt, af[c][0], af[c][1], tpos[0], tpos[-1], round(bits, 1), round(ident * 100, 1), hl, round(hi_, 1), hb, hnull, ins))

# consensus rows in the P26 frame: P26 itself, r6 and r5_27seqs where they align to it
def frame(nm):
    al = markx10(nm, cons[nm], P26, P); cell = ["-"] * L
    for q, t in pairs(al[2]):
        if t: cell[t - 1] = cons[nm][q - 1]
    return "".join(cell)
crows = [("CONS_r5_r6_P26 (r5-type head 2-137, r6 138-376)", P), ("CONS_r6_210seqs_on_P26", frame(R6)), ("CONS_r5_27seqs_on_P26", frame("r5_27seqs"))]
rows.sort(key=lambda r: ("ABCD".index(r[0]), r[1]))
WL = max([0] + [len(u) - (f - 1) for _, _, _, u, _, _, f, _ in rows]); WR = max([0] + [len(d) - (L - t) for _, _, _, _, _, d, _, t in rows])
with open(os.path.join(out, "r6_P26_review.aln.fa"), "w") as fo:
    for n, m in crows:
        fo.write(">%s\n%s%s%s\n" % (n.replace(" ", "_"), "-" * WL, m, "-" * WR))
    for g, _, n, u, m, d, f, t in rows:
        full = list("-" * WL + m + "-" * WR)
        for k, ch in enumerate(reversed(u)): full[WL + f - 2 - k] = ch            # ends at the column before P26 position f
        for k, ch in enumerate(d): full[WL + t + k] = ch                          # starts at the column after P26 position t
        fo.write(">%s\n%s\n" % (n, "".join(full)))
with open(os.path.join(out, "r6_P26_review.tsv"), "w") as fo:
    fo.write("group\tcopy\tvotes\tbest\tstatus\tP26_from\tP26_to\tbits_vs_P26\tident_vs_P26\thead_bp\thead_ident\thead_bits\thead_null_bits\tins\n")
    for r in info: fo.write("\t".join(map(str, r)) + "\n")
with open(os.path.join(out, "summary.tsv"), "w") as fo:
    fo.write("group\tcopies\tflanks_unique\tfs9_status\tshown\tP26_from_median\thead_above_null\thead_ge_50bp_above_null\n")
    for g in "ABCD":
        sel = [r for r in info if r[0] == g]
        fr = sorted(r[5] for r in sel)
        above = sum(1 for r in sel if r[11] > r[12]); big = sum(1 for r in sel if r[11] > r[12] and r[9] >= 50)
        fo.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\n" % (g, counts[g][0], counts[g][1], dict(counts[g][2]), len(sel), fr[len(fr) // 2] if fr else "-", above, big))
print(open(os.path.join(out, "summary.tsv")).read())
