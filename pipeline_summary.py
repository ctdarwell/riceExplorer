import pandas as pd, glob, numpy as np, sys

files = glob.glob1(sys.argv[1], '*_assocSNPs.csv') #assocSNPs is all extreme haplos
in_grp = sys.argv[2].replace('_', ' ').split(',')

hi, lo = pd.DataFrame(), pd.DataFrame()

for f in files:
    if '_lo_' in f: lo = pd.concat([lo, pd.read_csv(f"{sys.argv[1]}/{f}", header=0)], axis=0)
    if '_hi_' in f: hi = pd.concat([hi, pd.read_csv(f"{sys.argv[1]}/{f}", header=0)], axis=0)

try:
    hi_nfv = hi[hi.type.isin(in_grp)]
    hi_alt = hi.locus[~hi.type.isin(in_grp)].unique().tolist()
except:
    pass

try:
    lo_nfv = lo[lo.type.isin(in_grp)]
    lo_alt = lo.locus[~lo.type.isin(in_grp)].unique().tolist()
except:
    pass

hi_nfv[~hi_nfv.locus.isin(hi_alt)].drop_duplicates().to_csv(f'{sys.argv[1]}/hiVal_uniqueResults.csv', index=False)
lo_nfv[~lo_nfv.locus.isin(lo_alt)].drop_duplicates().to_csv(f'{sys.argv[1]}/loVal_uniqueResults.csv', index=False)


if not hi.empty:
    hi = hi.drop_duplicates()
    hi.to_csv(f"{sys.argv[1]}/fullTable_highValGenos.csv", index=False)
    hi = hi.acc.value_counts().reset_index().rename(
           columns={'index': 'acc', 'acc': 'count'})
    tmp = pd.DataFrame(np.array(['hi'] * hi.shape[0]).reshape(hi.shape[0],1), columns = list("X"))
    hi = pd.concat([hi, tmp], axis = 1)

if not lo.empty:
    lo = lo.drop_duplicates()
    lo.to_csv(f"{sys.argv[1]}/fullTable_lowValGenos.csv", index=False)
    lo = lo.acc.value_counts().reset_index().rename(
           columns={'index': 'acc', 'acc': 'count'})
    tmp = pd.DataFrame(np.array(['lo'] * lo.shape[0]).reshape(lo.shape[0],1), columns = list("X"))
    lo = pd.concat([lo, tmp], axis = 1)

df = pd.concat([hi, lo], axis = 0)

df.columns = ['accession','count','phenotype']

df.to_csv(f"{sys.argv[1]}/ALL_accession_cnts.csv", index=False)
