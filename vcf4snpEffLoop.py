import pandas as pd, sys, glob, os

vcfs = sys.argv[1] #e.g. LOC*vcf or rgd_*vcf
fldr = sys.argv[2]
panel = sys.argv[3] #eg IRIS W00
files = glob.glob1(fldr, vcfs)
cnt = 0

for vcf in files:
    cols = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','tbr']

    samps, dat = [], []
    for line in open(f"{fldr}/{vcf}"):
        if line.startswith(panel): samp = line.strip()
        else:
            try:
                samps.append(samp)
            except:
                continue
            dat.append(line.strip().split('\t'))

    if len(dat) == 0: #THIS MAYBE INEFFECTIVE
        try:
            del cols, line, samp
            os.remove(f"{fldr}/{vcf}")
        except: continue
        cnt += 1
        continue

    df = pd.DataFrame(dat)
    df = pd.concat([df, pd.DataFrame(samps)], axis=1)
    cols.append('samp')
    df.columns = cols
    del df['FORMAT']
    del df['tbr']

    df.INFO = '.'
    df.QUAL = '.'
    df.FILTER = '.'
    df.POS = df.POS.astype(int)
    df = df.sort_values('POS')

    #write csv of posn samples
    posn_samps = []
    for pos in df.POS.unique(): posn_samps.append([pos, df.samp[df.POS == pos].tolist()])
    posn_samps = pd.DataFrame(posn_samps)
    posn_samps.columns = ['marker', 'samples']
    posn_samps.to_csv(f"{fldr}/{vcf.replace('.vcf','.posn_samps')}.csv", index = False)

    del df['samp']
    df = df.drop_duplicates() #MAY NEED TO THINK ABOUT THIS!!!

    df.to_csv(f"{fldr}/{vcf}", index = False, sep='\t')
    del cols, line, samp

#The following is reqd at L30 if working with half processed folder

    '''
    if 'tbr' not in cols: #count no. of files already modifies (debugging)
        try:
            del cols, line, samp
            os.remove(f"{fldr}/{vcf}")
            print("BLX")
        except: continue
        cnt += 1
        continue
    '''
