import pandas as pd, sys

file = sys.argv[1] #main sample file
dat = pd.read_csv(file, header=0) #all IRRI samps w this trait

try: #remove any countries not wanted to be considered [optional]
    excl = sys.argv[2].replace('_', ' ').split(',')
    dat = dat[~dat.country.isin(excl)]
except: pass

try: #remove any individual samples not wanted to be considered [optional]
    tbr = pd.read_csv(sys.argv[3], header=0)
    dat = dat[~dat.acc.isin(tbr.acc.tolist())]
except: pass

df = pd.DataFrame()
for reg in dat.region.unique():
    tmp = dat[dat.region == reg]
    for sub in tmp.type.unique():
        tab = tmp[tmp.type == sub].sort_values('trait')
        df = pd.concat([df, tab.iloc[0, :]], axis = 1)
        if tab.at[tab.index[0], 'trait'] != tab.at[tab.index[-1], 'trait'] > 1: df = pd.concat([df, tab.iloc[-1, :]], axis = 1)

df = df.T #reqd from all DB samples with phenotypes

df.to_csv('samples_to_download.csv', index=False)

