import pandas as pd

trait = "GL"
dat = pd.read_csv('irri_grain_length.csv', header=0) #all IRRI samps w this trait
ttf = pd.read_csv('new_samps.csv', header=None)[0].tolist() #NB this contains list of all samps in new_irri_samps fldr

#CHECK COLUMN HEADERS - i.e. VAR, REGION & trait name
try: del dat['NAME'], dat['ACCESSION']
except: pass
dat.columns = ['IRIS_ID','SUBPOPULATION', 'VAR', 'COUNTRY','REGION', trait]
dat = dat[dat.COUNTRY != 'Thailand']

irri72 = pd.read_csv('irri_72vars_tidied.csv', header=0)
irri = irri72.DNA_ID[irri72.Country_Origin != 'Thailand'].tolist()
irri_samps = dat.loc[dat.IRIS_ID.isin(irri)] #gives non-Thai core samps with trait data

df = pd.DataFrame()
for reg in dat.REGION.unique():
    tmp = dat[dat.REGION == reg]
    for sub in tmp.VAR.unique():
        tab = tmp[tmp.VAR == sub].sort_values(trait)
        df = pd.concat([df, tab.iloc[0, :]], axis = 1)
        if tab.at[tab.index[0], trait] != tab.at[tab.index[-1], trait] > 1: df = pd.concat([df, tab.iloc[-1, :]], axis = 1)

df = df.T #reqd from all DB samples with phenotypes

#remove samples already downloaded from irri
df = df.loc[~df.IRIS_ID.isin(irri)] #combined w irri_samps gives all samps for anal
print(df)

#remove samples already downloaded from prev anals
df2 = df.loc[~df.IRIS_ID.isin(ttf)]
print(df2)

#files needed to be uploaded
#df2.to_csv(f"irri_{trait}_toAdd.csv", index=False)

#mk paths & data file from irri_samps + df
tmp = pd.concat([irri_samps, df])
#tmp.to_csv(f"irri{tmp.shape[0]}_{trait}_data.csv", index=False)
tmp.GL.hist()
'''
g = open(f"irri{tmp.shape[0]}_{trait}_paths.csv", 'w')
for s in irri_samps.IRIS_ID: g.write(f"../orig_irri_nipp/{s}.rechrom.vcf.gz\n")
for s in df.IRIS_ID: g.write(f"../new_irri_nipp/{s}.rechrom.vcf.gz\n")
g.close()
'''
