#!/bin/bash
# Rebuild of the Tal sicista_mito.html alignments on the read-POLISHED Sb1 mtDNA (2026-10-06).
# Same steps as aln/ (ali1.sh, ali2.sh): only the Sb1 backbone changes. ptg633_unit2 (the assembly unit) has
# ONT indel errors that frameshift COX1, COX3, ND5 (read-proven) and ND2 (probable); see pseudo_check/.
# The S. trizona NUMT alignment (backbone OZ418355.1) does not involve the Sb1 sequence and is not rebuilt.
set -euo pipefail
export PATH=$HOME/Sicista2026/envs/sicista_asm/bin:$HOME/miniforge3/envs/sinederella/bin:$PATH TMPDIR=$HOME/tmp
M=$HOME/Sicista2026/mito
mkdir -p $M/aln_v2
cd $M/aln_v2 || exit 1
awk 'NR==1{print ">Sb1_own_mtDNA_polished"; next} {print}' ../pseudo_check/seqs/unit2.fa > polished.fa
# 1. mitogenomes: the old all_mito.fa with the Sb1 row replaced
python3 - <<'PY'
rows=[];n=None
for l in open('../aln/all_mito.fa'):
    l=l.rstrip('\n')
    if l.startswith('>'): n=[l[1:],[]]; rows.append(n)
    else: n[1].append(l)
pol=''.join(l.strip() for l in open('polished.fa') if not l.startswith('>'))
out=open('all_mito.fa','w'); k=0
for h,s in rows:
    if 'Sb1_own_mtDNA' in h: h='Sbet_Sb1_own_mtDNA_polished'; s=[pol]; k+=1
    out.write('>%s\n%s\n'%(h,''.join(s)))
assert k==1, k
print(len(rows),'mitogenomes')
PY
mafft --quiet --thread 8 --auto --adjustdirection all_mito.fa > sicista_mitogenomes.aln.fa
# 2. control region (after tRNA-Pro, NC_069019 pos 15465), realigned alone
python3 - <<'PY'
def rd(f):
    out=[];n=None
    for l in open(f):
        l=l.rstrip('\n')
        if l.startswith('>'): n=[l[1:],[]]; out.append(n)
        else: n[1].append(l)
    return [(h,''.join(s)) for h,s in out]
A=rd('sicista_mitogenomes.aln.fa')
ref=[s for h,s in A if 'NC_069019' in h][0]
cols=[i for i,c in enumerate(ref) if c!='-']
c0=cols[15464]
with open('cr_raw.fa','w') as o:
    for h,s in A: o.write('>%s\n%s\n'%(h.replace('_R_',''),s[c0:].replace('-','')))
PY
mafft --quiet --thread 8 --localpair --maxiterate 100 --adjustdirection cr_raw.fa > sicista_control_region.aln.fa
# 3. Sb1 mt-like contig segments and nuclear copies on the polished backbone (fragments as before)
mafft --quiet --thread 8 --adjustdirection --keeplength --addfragments ../aln/mtlike_frag.fa polished.fa > sb1_mtlike_contigs.aln.fa
mafft --quiet --thread 8 --adjustdirection --keeplength --addfragments ../aln/sb1_fam_frag.fa polished.fa > sb1_nuclear_lineageB_family.aln.fa
mafft --quiet --thread 8 --adjustdirection --keeplength --addfragments ../aln/sb1_other_frag.fa polished.fa > sb1_nuclear_other_numts.aln.fa
# 4. gene track on the polished sequence: GenBank NC_069019 features lifted through the new mitogenome alignment
python3 - <<'PY'
# lift_bed.py matches the first word of headers; MAFFT may prefix reverse-complemented rows with _R_ (neither row here)
PY
python3 lift_bed.py sicista_mitogenomes.NC_069019_genes.bed sicista_mitogenomes.aln.fa Sbet_NC_069019_DK_BMNT Sbet_Sb1_own_mtDNA_polished Sb1_own_mtDNA_polished > sb1_polished_genes.bed
sed 's/^Sb1_own_mtDNA_polished\t/Sbet_Sb1_own_mtDNA_polished\t/' sb1_polished_genes.bed > sb1_polished_genes.mitogenomes_row.bed
wc -l sb1_polished_genes.bed
seqkit stats *.aln.fa | cut -c1-100
grep -c '_R_' *.aln.fa || true
echo ALN_V2_DONE
