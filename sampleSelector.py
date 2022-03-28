import pandas as pd, sys

file = 'irri_grain_length.csv' #sys.argv[1] #main sample file
dat = pd.read_csv(file, header=0) #all IRRI samps w this trait

alt_cols = set(dat.columns.tolist()) - set(['acc', 'type', 'region', 'trait'])
samps_reqd = 200 #sys.argv[2]

try: #remove any countries not wanted to be considered [optional]
    excl = ['Thailand'] #sys.argv[3].replace('_', ' ').split(',')
    dat = dat[~dat.country.isin(excl)]
except: pass

try: #remove any individual samples not wanted to be considered [optional]
    tbr = pd.read_csv(sys.argv[4], header=0)
    dat = dat[~dat.acc.isin(tbr.acc.tolist())]
except: pass

#first take best and worst from region & accession type
df, done = [], []
for reg in dat.region.unique():
    tmp = dat[dat.region == reg]
    for sub in tmp.type.unique():
        tab = tmp[tmp.type == sub].sort_values('trait')
        df.append(tab.iloc[0, :].tolist()), done.append(tab.at[tab.index[0], 'acc'])
        if tab.at[tab.index[0], 'acc'] != tab.at[tab.index[-1], 'acc']:
            df.append(tab.iloc[-1, :].tolist()), done.append(tab.at[tab.index[-1], 'acc'])

done = []
#take best and worse values from other columns
if len(df) < samps_reqd:
    dat = dat[~dat.acc.isin(done)]
    dat = dat.sample(frac=1) #scramble df
    for col in alt_cols:
        done = []
        for var in dat[col].unique():
            tab = dat[dat[col] == var].sort_values('trait')
            df.append(tab.iloc[0, :].tolist()), done.append(tab.at[tab.index[0], 'acc'])
            if tab.at[tab.index[0], 'trait'] != tab.at[tab.index[-1], 'trait']:
                df.append(tab.iloc[-1, :].tolist()), done.append(tab.at[tab.index[-1], 'acc'])
            if len(df) >= samps_reqd: break

        dat = dat[~dat.acc.isin(done)]
        if len(df) >= samps_reqd: break


#if not enough take random sample
if  len(df) < samps_reqd:
    print('Finishing off with random sampling')
    dat = dat.sample(samps_reqd - len(df))
    for i in dat.index:
        df.append(dat.loc[i].tolist())

df = pd.DataFrame(df)
df.columns = dat.columns.tolist()
df.to_csv('samples_to_download.csv', index=False)

print(f"Identified {df.shape[0]} samples")

