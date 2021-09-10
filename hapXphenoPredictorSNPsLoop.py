import pandas as pd, sys, glob
import numpy as np
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

depvar = sys.argv[1]
data = sys.argv[2]
fldr = sys.argv[3]
files = glob.glob1(fldr, '*best')
dat = pd.read_csv(data, header=0)
dat = pd.concat([dat.IRIS_ID, dat[depvar]], axis=1) #core of df (which gets modified) must comprise only two columns

#alt CI methods: #1 gives for a single draw, #2 gives for the mean of N draws - https://stackoverflow.com/questions/15033511/compute-a-confidence-interval-from-sample-data
m, s = dat[depvar].mean(), dat[depvar].std()
#ci = stats.norm.interval(0.95, loc=m, scale=s)
ci = stats.norm.interval(0.95, loc=m, scale=s/np.sqrt(dat[depvar].shape[0]))

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
        dat[mk][dat.IRIS_ID.isin(samps)] = snps[4][snps[1] == mk].values[0] #alt alleles
        dat[mk][~dat.IRIS_ID.isin(samps)] = snps[3][snps[1] == mk].values[0] #ref alleles

    dat['hap'] = dat[mks].agg(''.join, axis=1) #aggregate haplotype columns into single column
    haps = dat.hap.value_counts().keys().tolist() #list haplotypes (is ordered by freq)
    dat['val'] = None
    for h in haps: dat.val[dat.hap == h] = dat[depvar][dat.hap == h].mean() #column of means by haplotype

    #id haplos that have mean trait value above or below 95% CI
    sign = pd.DataFrame()
    sign = dat[dat.groupby('hap')['hap'].transform('size') > 5][dat.val < ci[0]]
    sign = pd.concat([sign, dat[dat.groupby('hap')['hap'].transform('size') > 5][dat.val > ci[1]]])
    sign.val[sign.val < dat[depvar].mean()] = 'lo'
    sign.val[sign.val != 'lo'] = 'hi'
    dat = dat.iloc[:, :2]

    #tidy up df
    del sign[depvar], sign['hap']
    if sign.val.unique().shape[0] > 1: #if hi v lo comparion, rm any invariants - otherwise, leave all
        for mk in mks: #rm invariant SNPs in sign loci
            try:
                if sign[mk].unique().shape[0] == 1: del sign[mk]
            except: continue

    #write file if any interesting haplotypes left
    if not sign.empty: sign.to_csv(f"{fldr}/{msu.replace('ann.best','haps.csv')}", index = False)

