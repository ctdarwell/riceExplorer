import pandas as pd, sys, glob
import warnings
warnings.filterwarnings("ignore")

fldr = sys.argv[1]
targets = ['lo', 'hi']
files = glob.glob1(fldr, "*haps.csv")
ACCS = sys.argv[2]
accs = pd.read_csv(ACCS, header=0) 
types = sys.argv[3].replace('_', ' ').split(',')

for msu in files:
    try:
        haps = pd.read_csv(f"{fldr}/{msu}", header=0) #haps from IRRI samps w/ extreme values
        rgd = pd.read_csv(f"{fldr}/rgd_{msu.replace('haps.csv', 'vcf')}", header=0, sep='\t') #RGD samps of polymorphic sites
        samps_list = pd.read_csv(f"{fldr}/rgd_{msu.replace('haps.csv', 'posn_samps.csv')}", header=0) #list of sites with RGD samples possessing ALT alleles
        vcf = pd.read_csv(f"{fldr}/{msu.replace('haps.csv', 'vcf')}", header=0, sep='\t') #IRRI samps of polymorphic sites
    except: continue

    for target in targets:
        cnt, cnt2 = 1, 1
        haplos = haps.drop_duplicates(haps.columns[1:-1].tolist())[haps.val == target].iloc[:,1:-1] #unique IRRI haplotypes
        rgd = rgd[rgd.POS.isin(haplos.columns.tolist())] #Get rid of uninteresting RGD sites
        try: rgd['ALT'], rgd['MISC'] = rgd['ALT'].str.split(',', 1).str #clean the 'ALT' column (rm "<NON_REF>")
        except: pass

        novel, general = pd.DataFrame(), pd.DataFrame()
        for indx in haplos.index: #thru haplotypes of interest
            allor, samps = [], []
            for col in haplos.columns: #thru ea SNP in the haplotype
                try:
                    if haplos.loc[indx, col] == rgd.REF[rgd.POS == int(col)].values[0]: #Is it the reference allele at this site?
                        allor.append('REF')
                        tmp = samps_list.samples[samps_list.marker == int(col)].tolist()[0][2:-2].split("', '")
                        samps.append(list(set(accs.acc) - set(tmp))) #add all samples with REF allele at this site
                    else:
                        allor.append('ALT')
                        samps.append(samps_list.samples[samps_list.marker == int(col)].tolist()[0][2:-2].split("', '")) ##add all samples with ALT allele at this site
                except: #thrown if 'rgd' doesn't have polymorphic SNP from IRRI samps
                    if haplos.loc[indx, col] == vcf.REF[vcf.POS == int(col)].values[0]:
                        allor.append('REF') #as there is no polymorphism in RGD samples at this marker
                        samps.append(accs.acc.unique().tolist()) #add all samples as all have REF allele at this site
                    else: allor, samps = [], []
            try: int_samps = list(set.intersection(*map(set, samps))) #the intersection of the list of lists (samps) gives the indvs with that exact haplotype
            except: continue

            #check that ***ONLY*** NFVs of interest have the haplotype (and modify df)
            tab = accs[accs.acc.isin(int_samps)]
            if 'ALT' in allor and tab[~tab.type.isin(types)].empty: #id'd extreme variants UNIQUELY found in NFVs of interest
                tab['snps'] = '; '.join(haplos.columns.to_list())
                tab['alleles'] = '; '.join(haplos.loc[indx, :].to_list())
                tab['origin'] = '; '.join(allor)
                tab['locus'] = msu.replace('.haps.csv', '')
                tab['hap.no'] = cnt
                novel = pd.concat([novel, tab])
                cnt += 1

            #and catalog all extreme SNPs irrespective of variety type
            tab2 = accs[accs.acc.isin(int_samps)]
            if 'ALT' in allor: 
                tab2['snps'] = '; '.join(haplos.columns.to_list())
                tab2['alleles'] = '; '.join(haplos.loc[indx, :].to_list())
                tab2['origin'] = '; '.join(allor)
                tab2['locus'] = msu.replace('.haps.csv', '')
                tab2['hap.no'] = cnt2
                general = pd.concat([novel, tab2])
                cnt2 += 1

        if not novel.empty: novel.to_csv(f"{fldr}/{msu.replace('.haps.csv', '')}_potential_{target}_donors.csv", index=False)
        if not general.empty: general.to_csv(f"{fldr}/{msu.replace('.haps.csv', '')}_{target}_assocSNPs.csv", index=False)

