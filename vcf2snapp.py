import pandas as pd, sys, numpy as np

#For NJ calc etc consider: http://scikit-bio.org/docs/0.2.1/generated/skbio.tree.nj.html
#OR biopython: Phylogenetics with Bio.Phylo

vfile = sys.argv[1]
indir, outdir = sys.argv[2], sys.argv[3] 

#find start of VCF
cnt = 0
for line in open(f"{indir}/{vfile}"):
    if line.startswith('#CHROM'): break
    cnt += 1
    
vcf = pd.read_table(f"{indir}/{vfile}", header = cnt)
vcf.iloc[:, 9:] = vcf.iloc[:, 9:].applymap(lambda x: x.split(':')[0]) #use only allele info

#fnc to convert bin coded alleles to sequence codes
def myfunc(x): return dik1.get(x)
vfunc = np.vectorize(myfunc)

#convert bin coded alleles to sequence codes
arr = []
dik1 = {'./.': '?', '1/1': '2', '0/0': '0', '0/1': '1'}
for row in vcf.index:
    als = vfunc(np.array(vcf.iloc[row, 9:]))
    arr.append(als)

#extraxt sequences for each taxon from array
seqs = []
arr = np.array(arr)
for i in range(arr.shape[1]):
    seqs.append(''.join(arr[:, i]))

#write nexus seq data
def fas2nex(taxa, seqs):
    g = open(f"{outdir}/{vfile.replace(vfile.split('.')[-1], 'SNAPP.nex')}", "w")
    g.write(f"#NEXUS\nbegin data;\ndimensions ntax={len(seqs)} nchar={len(seqs[0])};\nFormat datatype=integerdata symbols='012' missing=?;\nmatrix\n")
    for tax in taxa:
        g.write(f"{tax}\n\t{seqs[taxa.index(tax)]}\n")
    g.write(";\nend;\n")
    g.close()


fas2nex(vcf.columns[9:].tolist(), seqs)



