import pandas as pd, sys, numpy as np

import glob

fldr = sys.argv[1]
vcfs = glob.glob1('.', '*faoCmbn.vcf') #sys.argv[1] #

#find start of VCF
cnt = 0
for line in open(f"{vcfs[0]}"):
    if line.startswith('#CHROM'): break
    cnt += 1
    
df = pd.read_table(f"{vcfs[0]}", header = cnt)
df.iloc[:, 9:] = df.iloc[:, 9:].applymap(lambda x: x.split(':')[0]) #use only allele info
df.iloc[:, 4:5] = df.iloc[:, 4:5].applymap(lambda x: x.replace(',<NON_REF>', '')) #use only allele info

for vcf in vcfs[1:]:
    #find start of VCF
    cnt = 0
    for line in open(vcf):
        if line.startswith('#CHROM'): break
        cnt += 1
        
    vcf = pd.read_table(vcf, header = cnt)
    vcf.iloc[:, 9:] = vcf.iloc[:, 9:].applymap(lambda x: x.split(':')[0]) #use only allele info
    vcf.iloc[:, 4:5] = vcf.iloc[:, 4:5].applymap(lambda x: x.replace(',<NON_REF>', '')) #use only allele info
    df = pd.concat([df, vcf])

df = df.reset_index(drop=True)

#rm Indels from REF AND ALT columns
lens = df.iloc[:, 4:5].applymap(lambda x: len(x)) #use only allele info
indxs = lens.index[lens.ALT > 1]
df = df.drop(indxs)
lens = df.iloc[:, 3:4].applymap(lambda x: len(x)) #use only allele info
indxs = lens.index[lens.REF > 1]
df = df.drop(indxs)


df.iloc[:, 9:] = df.iloc[:, 9:].applymap(lambda x: x.split(':')[0]) #use only allele info

#rm samples w >20% missing data
to_drop = []
x = df.shape[0]
for samp in df.columns[9:]:
    if df[samp].tolist().count('./.') / x > 0.2: 
        to_drop.append(samp)

df = df.drop(to_drop, axis=1)
df.to_csv('tmp.vcf', index=False, sep='\t')

g = open(f'{fldr}.allChroms.vcf', 'w')
g.write('##fileformat=VCFv4.1\n')
for line in open('tmp.vcf'):
    g.write(line)
g.close()
