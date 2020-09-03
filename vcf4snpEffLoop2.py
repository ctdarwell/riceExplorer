#reformat files outputted from BCFtools to snpEff readable
import pandas as pd, sys, glob, os

vcfs = sys.argv[1] #e.g. LOC*vcf or rgd_*vcf
fldr = sys.argv[2]
panel = sys.argv[3]
files = glob.glob1(fldr, vcfs)
cnt = 0

for vcf in files:
    samps, dat = [], []
    for line in open(f"{fldr}/{vcf}"): #create DF from BCFtools outputs
        if line.startswith(panel): samp = line.strip()
        else:
            try:
                samps.append(samp)
            except:
                cols = line.strip().split('\t')
            dat.append(line.strip().split('\t'))

    if len(dat) == 1:
        print(f"{vcf} - NO DATA - vcf4snpeff.py aborted!")
        try:
            del cols, line, samp
            os.remove(f"{fldr}/{vcf}")
        except: continue
        cnt += 1
        continue

    if 'misc' not in cols: #count no. of files already modifies (debugging)
        try:
            del cols, line, samp
            os.remove(f"{fldr}/{vcf}")
        except: continue
        cnt += 1
        continue

    #format DF to snpEff readable
    df = pd.DataFrame(dat[1:])
    df = pd.concat([df, pd.DataFrame(samps)], axis=1)
    cols.append('samp')
    df.columns = cols
    del df['FORMAT']
    del df['misc']

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
    df = df.drop_duplicates()

    df.to_csv(f"{fldr}/{vcf}", index = False, sep='\t')
    del cols, line, samp

print(f"vcf4snpEff - no. files already modified = {cnt}")

