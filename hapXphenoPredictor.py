import pandas as pd, sys, glob, numpy as np
import warnings
from scipy import stats
warnings.filterwarnings("ignore")

import concurrent.futures, itertools
max_workers = None #nProcesses - leave as default on HPC [ie None]; set as no. processors on your machine -1 for laptops/Desktops 
chunksize = 1 #break up of parallelisation chunks; large values may be faster for larger datasets

data = sys.argv[1] # 
fldr = sys.argv[2]
thresh = sys.argv[3]  # min/max vals eg 7,10 or 0.866 = 1.5sd, see https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule

try: #determine extreme phenotypic threshold (if expressed as abs vals)
    q1, q2 = float(thresh.split(',')[0]), float(thresh.split(',')[1])
    swch = 0
except:
    thresh = float(thresh)
    swch = 1

files = glob.glob1(fldr, '*best') #get snpEff output files
dat = pd.read_csv(data, header=0) #get sample name/phenotype data
dat = pd.concat([dat.acc, dat.trait], axis=1) #core of df (which gets modified) must comprise only two columns

#main processing function
def fnc(msu, dat):
    try: #load files sequentially (i.e. not parallelised but not time consuming)
        snps = pd.read_csv(f"{fldr}/{msu}", header = None, sep='\t')
        posn_samps = pd.read_csv(f"{fldr}/{msu.replace('ann.best','posn_samps.csv')}", header=0)
    except: return

    #build haplotypes at separate columns in df
    mks = np.sort(snps[1]) #NB no headers with snpEff outputs
    for mk in mks:
        dat[mk] = np.nan
        tmp = posn_samps[posn_samps.marker == mk]
        samps = tmp.samples.tolist()[0][2:-2].split("', '") #bit ugly but works from outputs
        dat[mk][dat.acc.isin(samps)] = snps[4][snps[1] == mk].values[0] #alt alleles
        dat[mk][~dat.acc.isin(samps)] = snps[3][snps[1] == mk].values[0] #ref alleles

    dat['hap'] = dat[mks].agg(''.join, axis=1) #aggregate haplotype columns into single column
    haps = dat.hap.value_counts().keys().tolist() #list haplotypes (is ordered by freq)
    dat['val'] = None
    for h in haps: dat.val[dat.hap == h] = dat.trait[dat.hap == h].mean() #column of means by haplotype

    #determine extreme phenotypic threshold (if expressed as a ci percentile)
    if swch == 1:
        m, s = dat.val.unique().mean(), dat.val.unique().std()
        ci = stats.norm.interval(thresh, loc=m, scale=s)
        q1, q2 = ci[0], ci[1]

    #id haplos that have mean trait value above or below thresh
    sign = pd.DataFrame()
    sign = dat[dat.groupby('hap')['hap'].transform('size') > 5][dat.val < q1] #lower percentiles
    sign = pd.concat([sign, dat[dat.groupby('hap')['hap'].transform('size') > 5][dat.val > q2]]) #upper percentiles
    sign.val[sign.val < dat.trait.mean()] = 'lo'
    sign.val[sign.val != 'lo'] = 'hi'
    dat = dat.iloc[:, :2] #strip dat back to original two columns for next input file

    #tidy up df
    del sign['trait'], sign['hap']
    if sign.val.unique().size > 1: #if hi v lo comparion, rm any invariants - otherwise, leave all
        for mk in mks: #rm invariant SNPs in sign loci
            try:
                if sign[mk].unique().size == 1: del sign[mk]
            except: continue

    #write file if any interesting haplotypes left
    if not sign.empty: sign.to_csv(f"{fldr}/{msu.replace('ann.best','haps.csv')}", index = False)

def main():
    #run parallelised function
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        executor.map(fnc, files, itertools.repeat(dat), chunksize=chunksize)

if __name__ == '__main__': main()

