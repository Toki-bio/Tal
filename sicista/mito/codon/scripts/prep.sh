#!/bin/bash
# Inputs for the Sicista mitochondrial pseudogene check (2026-10-06)
set -euo pipefail
cd ~/Sicista2026/mito/pseudo_check || exit 1
export PATH=$HOME/miniforge3/envs/sinederella/bin:$PATH TMPDIR=$HOME/tmp
mkdir -p seqs gb
ASM=$HOME/Sicista2026/assembly/primary_eval/hifiasm_primary.p_ctg.fa
TRI=$HOME/Sicista2026/reference/GCA_982266845.1_mSicTri.1_genomic.fna
# 1. GenBank records (61), renamed to the labels used on the Tal page; full flat files for the annotation check
python3 - <<'PY'
import re
lab={}
for l in open('labels.txt'):
    l=l.strip(); m=re.search(r'(NC_\d{6}|[A-Z]{1,2}\d{5,})',l)
    if m: lab[m.group(1)]=l
out=open('seqs/gb.fa','w'); ids=[]
for line in open('../genbank/all_sicista_mito.fa'):
    if line.startswith('>'):
        acc=line[1:].split()[0]; ids.append(acc); base=acc.split('.')[0]
        out.write('>'+lab.get(base,'GB_'+base)+'\n')
    else: out.write(line)
open('gb/ids.txt','w').write('\n'.join(ids)+'\n')
PY
grep -c '>' seqs/gb.fa; grep '>GB_' seqs/gb.fa || echo all_labelled
[ -s gb/all.gb ] || curl -sS "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=gb&retmode=text&id=$(paste -sd, gb/ids.txt)" > gb/all.gb
grep -c '^LOCUS' gb/all.gb
# 2. own mtDNA (lineage A, read-verified)
# read-polished (samtools consensus -X r10.4_sup on reads/A_unit2.bam): fixes 3 ONT frameshift errors of the assembly (COX1 +CT, COX3 +TTTTC, ND5 -ACTC) and +AAA at 5192
# + ND2 C-run at 4896: reads C5 53% / C6 32% / C4 15% (ONT homopolymer undercall); C6 in the 33 GenBank records with the same flanks (MZ570947 has C5) and
#   C5 frameshifts ND2 -> C6 used (PROBABLE, not read-proven; needs Sanger/Illumina)
python3 - <<'PY'
s=''.join(l.strip() for l in open('unit2_consensus.fa') if not l.startswith('>'))
assert s[4885:4900]=='AACATTTTAACCCCC' and s[4900]!='C', s[4885:4905]
s=s[:4895]+'C'+s[4895:]
open('seqs/unit2.fa','w').write('>Sb1_own_mtDNA_polished\n'+'\n'.join(s[i:i+60] for i in range(0,len(s),60))+'\n')
PY
# gene coordinates lifted from units/genes_unit2.tsv by the five indels (assembly coords; 4895 C-run, 5192 +AAA, 6639 +CT, 8895 +TTTTC, 12471 -ACTC)
awk -F'	' 'BEGIN{OFS="	"} {s=$2; e=$3; s+=(s>4895?1:0)+(s>5192?3:0)+(s>6639?2:0)+(s>8895?5:0)-(s>12474?4:0); e+=(e>4895?1:0)+(e>5192?3:0)+(e>6639?2:0)+(e>8895?5:0)-(e>12474?4:0); print $1,s,e}' ../units/genes_unit2.tsv > genes_polished.tsv
# 3. betulina nuclear copies: 82 lineage-B family (span>=14 kb, ~92% to own mtDNA) + young own-lineage copies
# loci padded by 1.5 kb each side, so gene pieces at the insertion junctions are not cut off by the locus coordinates
awk -F'\t' 'NR>1 && $5>=14000 {s=$2-1500; if(s<1)s=1; print $1":"s"-"$3+1500"\t"$4"\t"int($7)}' ../bet_numt/bet_numt_loci.tsv > seqs/bet_loci.txt
printf 'ptg000880l:1-13509\tplus\t99\n' >> seqs/bet_loci.txt
: > seqs/bet_numt.fa
while IFS=$'\t' read -r reg str pid; do
  tag=NUMTB; [ "$pid" -ge 98 ] && tag=NUMTA; [ "$pid" -lt 90 ] && tag=NUMTold
  name="Sb1_${tag}_${reg//[:-]/_}"
  if [ "$str" = minus ]; then samtools faidx -i "$ASM" "$reg"; else samtools faidx "$ASM" "$reg"; fi | awk -v n="$name" 'NR==1{print ">"n; next} {print}' >> seqs/bet_numt.fa
done < seqs/bet_loci.txt
grep '>' seqs/bet_numt.fa | cut -d_ -f2 | sort | uniq -c
# 4. S. trizona nuclear loci >= 1 kb
awk -F'\t' 'NR>1 && $5>=1000 {s=$2-1500; if(s<1)s=1; print $1":"s"-"$3+1500"\t"$4"\t"int($7)}' ../trizona_numt/tri_numt_loci.tsv > seqs/tri_loci.txt
: > seqs/tri_numt.fa
while IFS=$'\t' read -r reg str pid; do
  name="Stri_NUMT${pid}_${reg//[:.-]/_}"
  if [ "$str" = minus ]; then samtools faidx -i "$TRI" "$reg"; else samtools faidx "$TRI" "$reg"; fi | awk -v n="$name" 'NR==1{print ">"n; next} {print}' >> seqs/tri_numt.fa
done < seqs/tri_loci.txt
grep -c '>' seqs/tri_numt.fa
cat seqs/unit2.fa seqs/gb.fa seqs/bet_numt.fa seqs/tri_numt.fa > seqs/all.fa
seqkit stats seqs/*.fa
echo PREP_DONE
