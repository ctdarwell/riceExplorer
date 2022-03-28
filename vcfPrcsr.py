import pandas as pd, sys

chrom = sys.argv[1] # 

vfile = f"Chr{chrom}.Mge.vcf"
max_missing = 0.2

val = int((1 - max_missing) * 100)
snps = pd.read_csv(f'Chr{chrom}.neworder.txt', header=None)

snps.columns = ['snps']
snps.iloc[:, :1] = snps.iloc[:, :1].applymap(lambda x: x.split(':')[1]) #use only allele info
locs = snps.snps.astype(int).tolist()


#find start of VCF
cnt = 0
for line in open(vfile):
    if line.startswith('#CHROM'): break
    cnt += 1

vcf = pd.read_table(vfile, header = cnt, sep='\t')

#vcf = vcf[~vcf.FILTER.str.startswith('LowQual')]
vcf = vcf.reset_index(drop=True)

#rm SNP with >1 record
cnts = vcf.groupby(['POS']).size().reset_index().rename(columns={0:'cnt'})
rm = cnts.POS[cnts.cnt > 1].tolist()
vcf = vcf[~vcf.POS.isin(rm)]

#keep only from id'd SNPs
vcf[vcf.POS.isin(locs)].to_csv('tmp.vcf', index=False, sep='\t')

#rm multi-allelic loci
g = open('tmp2.vcf', 'w')
for line in open('tmp.vcf'):
    if '/2' not in line and '/3' not in line and '/4' not in line and '/5' not in line and '/6' not in line and '/7' not in line and '/8' not in line and '/9' not in line and '/10' not in line:
        g.write(line)
g.close()


#load edited VCF
cnt = 0
for line in open('tmp2.vcf'):
    if line.startswith('#CHROM'): break
    cnt += 1
    
vcf = pd.read_table('tmp2.vcf', header = cnt, sep='\t')

vcf2 = vcf.copy()
vcf2.iloc[:, 9:] = vcf2.iloc[:, 9:].applymap(lambda x: x.split(':')[0]) #use only allele info

#missing data <20%
reduc = vcf2[vcf2[vcf2 == './.'].count(axis=1) < (vcf2.shape[1] * max_missing)]
vcf = vcf.iloc[reduc.index]

vcf.to_csv(vfile.replace('Mge', 'faoCmbn'), index=False, sep='\t')
vcf2 = vcf.copy()
vcf2.iloc[:, 9:] = vcf2.iloc[:, 9:].applymap(lambda x: x.split(':')[0]) #use only allele info

#rm samples w >20% missing data
to_drop = []
x = vcf2.shape[0]
for samp in vcf2.columns[9:]:
    if vcf2[samp].tolist().count('./.') / x > 0.2: 
        to_drop.append(samp)

faoLD = vcf.drop(to_drop, axis=1)

faoLD.to_csv('tmp.vcf', index=False, sep='\t')

g = open(vfile.replace('Mge', 'faoLD'), 'w')
g.write('##fileformat=VCFv4.1\n')
for line in open('tmp.vcf'):
    g.write(line)
g.close()

