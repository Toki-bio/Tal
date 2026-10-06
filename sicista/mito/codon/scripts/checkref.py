# translate the reference genes (vertebrate mito code) - all must be stop-free inside
b="TCAG"; aa="FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG"; code={}; k=0
for x in b:
    for y in b:
        for z in b: code[x+y+z]=aa[k]; k+=1
code.update({"TGA":"W","ATA":"M","AGA":"*","AGG":"*"})
bad=0
for g in ["ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB"]:
    s=[l.strip() for l in open(f"genes/{g}.query.fa")][1]
    p="".join(code.get(s[i:i+3],"X") for i in range(0,len(s)-2,3))
    n=p[:-1].count("*"); bad+=n
    print(g,len(s),"internal stops",n,"start",s[:3],"last codon/partial",s[len(s)-len(s)%3-3:] if len(s)%3 else s[-3:])
print("REF_OK" if bad==0 else "REF_HAS_STOPS")
