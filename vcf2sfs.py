import pandas as pd, sys, numpy as np
import matplotlib.pyplot as plt
#ARE THERE UNK CHROMS??

vfile = sys.argv[1]
indir, outdir = sys.argv[2], sys.argv[3] 

#find start of VCF
cnt = 0
for line in open(f"{indir}/{vfile}"):
    if line.startswith('#CHROM'): break
    cnt += 1
    
vcf = pd.read_table(f"{indir}/{vfile}", header = cnt)

vcf['#CHROM'] = vcf['#CHROM'].astype(str)

vcf.iloc[:, 9:] = vcf.iloc[:, 9:].applymap(lambda x: x.split(':')[0]) #use only allele info

#convert SNP codes to summable floats
def myfunc(x): return dik1.get(x)
vfunc = np.vectorize(myfunc)
dik1 = {'./.': np.nan, '1/1': 2.0, '0/0': 0.0, '0/1': 1.0}
arr = vfunc(np.array(vcf.iloc[:, 9:]))
arr2 = np.nan_to_num(arr) #summing any NaNs along axis generates a nan value!

#calc freqs/site
al_freqs = np.apply_along_axis(sum, 1, arr2) / ((arr.shape[1] * 2) - np.apply_along_axis(sum, 1, np.isnan(arr)))

#plot SFSs
for chrom in vcf['#CHROM'].unique():
    indxs = vcf.index[vcf['#CHROM'] == chrom]
    plt.hist(al_freqs[indxs], bins=100, color='k')
    plt.xlabel("Reference allele frequency distribution")
    plt.ylabel("Freq.")
    plt.savefig(f"{outdir}/{vfile.replace('.vcf','')}_Chrom{chrom}_SFS.pdf", dpi = 450)

plt.hist(al_freqs, bins=100, color='k')
plt.xlabel("Site reference allele frequencies")
plt.ylabel("Freq.")
plt.savefig(f"{outdir}/{vfile.replace('.vcf','')}_allSites_SFS.pdf", dpi = 450)

tab = pd.concat([vcf[['#CHROM','POS']], pd.DataFrame(al_freqs)], axis=1)
tab.columns = ['#CHROM','POS','ref_frq']
tab.to_csv(f"{outdir}/{vfile.replace('.vcf','_SFStab.csv')}", index=False)


