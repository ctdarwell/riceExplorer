import pandas as pd, sys, glob
import numpy as np
import warnings
warnings.filterwarnings("ignore")

data = sys.argv[1] #'irri100_GL_data.csv' #
fldr = sys.argv[2]
files = glob.glob1(fldr, '*best')
dat = pd.read_csv(data, header=0)
dat = pd.concat([dat.acc, dat.trait], axis=1) #core of df (which gets modified) must comprise only two columns

q = float(sys.argv[3]) #percentile delimiting extreme phenotypes
q1 = dat.quantile(q=(1-q))[0]
q2 = dat.quantile(q=q)[0] 

for msu in files:
    try:
        snps = pd.read_csv(f"{fldr}/{msu}", header = None, sep='\t')
        posn_samps = pd.read_csv(f"{fldr}/{msu.replace('ann.best','posn_samps.csv')}", header=0)
    except:
        print(f"{msu} - has no potentially significant SNPs") #if no output was produced earlier
        continue

    #build haplotypes at separate columns in df
    mks = np.sort(snps[1])
    for mk in mks:
        dat[mk] = np.nan
        tmp = posn_samps[posn_samps.marker == mk]
        samps = tmp.samples.tolist()[0][2:-2].split("', '")
        dat[mk][dat.acc.isin(samps)] = snps[4][snps[1] == mk].values[0] #alt alleles
        dat[mk][~dat.acc.isin(samps)] = snps[3][snps[1] == mk].values[0] #ref alleles

    dat['hap'] = dat[mks].agg(''.join, axis=1) #aggregate haplotype columns into single column
    haps = dat.hap.value_counts().keys().tolist() #list haplotypes (is ordered by freq)
    dat['val'] = None
    for h in haps: dat.val[dat.hap == h] = dat.trait[dat.hap == h].mean() #column of means by haplotype

    #id haplos that have mean trait value above or below 95% CI
    sign = pd.DataFrame()
    sign = dat[dat.groupby('hap')['hap'].transform('size') > 5][dat.val < q1] #lower percentiles
    sign = pd.concat([sign, dat[dat.groupby('hap')['hap'].transform('size') > 5][dat.val > q2]]) #upper percentiles
    sign.val[sign.val < dat.trait.mean()] = 'lo'
    sign.val[sign.val != 'lo'] = 'hi'
    dat = dat.iloc[:, :2]

    #tidy up df
    del sign['trait'], sign['hap']
    if sign.val.unique().shape[0] > 1: #if hi v lo comparion, rm any invariants - otherwise, leave all
        for mk in mks: #rm invariant SNPs in sign loci
            try:
                if sign[mk].unique().shape[0] == 1: del sign[mk]
            except: continue

    #write file if any interesting haplotypes left
    if not sign.empty: sign.to_csv(f"{fldr}/{msu.replace('ann.best','haps.csv')}", index = False)

