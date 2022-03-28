import pandas as pd, sys, glob


search=sys.argv[1]
fldr = '.'
files = glob.glob1(fldr, '*faoLD.vcf')

pops = pd.read_csv('pops.csv', header=0)

for vfile in files:
    #find start of VCF
    cnt = 0
    for line in open(vfile):
        if line.startswith('#CHROM'): break
        cnt += 1
        
    vcf = pd.read_csv(vfile, header = cnt, sep='\t')
    
    vcf.iloc[:, 9:] = vcf.iloc[:, 9:].applymap(lambda x: x.split(':')[0]) #use only allele info
    
    tmp = pops.samps[pops.population.str.contains(search)].tolist()
    
    ind = sorted(list(set(vcf.columns.tolist()).intersection(tmp)))
    
    vcf[vcf.columns[:9].tolist() + ind].to_csv('tmp.vcf', index=False, sep='\t')
    
    g = open(vfile.replace('faoLD', f'{search}.LD'), 'w')
    g.write('##fileformat=VCFv4.1\n')
    for line in open('tmp.vcf'):
        g.write(line)
    g.close()



