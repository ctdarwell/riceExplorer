import pandas as pd, glob, numpy as np, sys, itertools

fldr = sys.argv[1] #$fldr
ld_fldr = sys.argv[2] #'.'
srch = sys.argv[3] #$gene_id
msu = sys.argv[4] #$genes
known = sys.argv[5] #$idd_genes

data = pd.read_csv(f'{fldr}/fullTable_highValGenos.csv', header=0)
nfvs = data[(data.type == 'Landrace') | (data.type == 'RD variety')]
ecs = data.locus[data.type == 'RGD improved line'].unique().tolist()
unq_nfvs = nfvs.locus[~nfvs.locus.isin(ecs)].drop_duplicates().tolist()

msu = pd.read_csv(msu, header=0)
try: gl = pd.read_csv(known, header=None)[0].tolist() # maybe not included
except: gl = []

all_locs = np.array(gl + unq_nfvs)

lds = glob.glob1(ld_fldr, '*IND.LD.csv')

for ld in lds:
    df = []
    tab = pd.read_csv(f'{ld_fldr}/{ld}', header=0)
    tab.index = tab["Unnamed: 0"]
    tab = tab.iloc[:, 1:]
    indxs = tab.columns.astype(int)
    locs = all_locs[np.where(np.char.startswith(all_locs, f"{srch}{ld[3:5]}"))]
    for it in itertools.combinations(locs, 2):
        l1 = [msu.start[msu.locus == it[0]].iloc[0], msu.stop[msu.locus == it[0]].iloc[0]]
        l2 = [msu.start[msu.locus == it[1]].iloc[0], msu.stop[msu.locus == it[1]].iloc[0]]
        pos1 = indxs[np.where((indxs > l1[0]) & (indxs < l1[1]))]
        pos2 = indxs[np.where((indxs > l2[0]) & (indxs < l2[1]))]
        if pos1.size > 0 and pos2.size > 0: #if LD has been recorded from a SNP within the gene range
            m, M = min(pos2[0], pos1[0]), max(pos2[0], pos1[0])
            df.append([it[0], it[1], tab.at[m, str(M)]])
        else: #if not - use the nearest SNP either side of the gene
            pos1a = indxs[np.where(indxs < l1[0])]
            pos1b = indxs[np.where(indxs > l1[1])]
            pos2a = indxs[np.where(indxs < l2[0])]
            pos2b = indxs[np.where(indxs > l2[1])]
                
            #check if both SNPs are unrepresented in LD matrix
            if pos1.size > 0: p1 = pos1[0]
            elif (l1[0] - pos1a[-1]) < (pos1b[0] - l1[1]): p1 = pos1a[-1]
            else: p1 = pos1b[0]
            if pos2.size > 0: p2 = pos2[0]
            elif (l2[0] - pos2a[-1]) < (pos2b[0] - l2[1]): p2 = pos2a[-1]
            else: p2 = pos2b[0]
            
            m, M = min(p2, p1), max(p2, p1)
            df.append([it[0], it[1], tab.at[m, str(M)]])

    df = pd.DataFrame(df)
    df.columns = ['Loc1','Loc2','LD']
    df.to_csv(ld.replace('csv','LGHtab.csv'), index=False)
