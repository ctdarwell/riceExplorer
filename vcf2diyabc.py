import pandas as pd, sys, numpy as np

vfile =  sys.argv[1]
indir, outdir = sys.argv[2], sys.argv[3]
pops = pd.read_csv(f"{sys.argv[2]}/{sys.argv[4]}", header=0) #samps, pops

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
dik1 = {'./.': '9', '1/1': '2', '0/0': '0', '0/1': '1'}
for row in vcf.index:
    als = vfunc(np.array(vcf.iloc[row, 9:]))
    arr.append(list(als))

arr = np.array(arr)
dik = dict(zip(pops.samps, pops.population))

#write nexus seq data
def vcf2diy(taxa, seqs):
    g = open(f"{outdir}/{vfile.replace(vfile.split('.')[-1], 'DIYABC.txt')}", "w")
    g.write(f"{vfile} to diyABC format <NM=1.0NF> <MAF=hudson>\nIND SEX POP")
    for i in range(seqs.shape[0]): g.write(" A")
    g.write('\n')
    for tax in taxa:
        g.write(f"{tax} 9 {dik.get(tax)} {' '.join(seqs[:, taxa.index(tax)])}\n")
    g.close()

vcf2diy(vcf.columns[9:].tolist(), arr)

