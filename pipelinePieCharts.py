import pandas as pd, sys, glob, cv2, numpy as np, os
import matplotlib.pyplot as plt

#NB THERE IS A WEIRD ISSUE RUNNING THIS LOOPED IN ANACONDA PROMPT WITH REGARD PIE SIZES: crops.append([pie.split('_')[0], cv2.resize(pic, (0,0), fx=0.8, fy=0.8)])
#IMPORTANT: FOR HI RES IMG, COMMENT OUT LAST LINE (for png in glob.glob1(".", f"*{target}_pie.png"): os.remove(png))
#THEN USE INDV PIES IN INKSCAPE, CONVERT TO PDF, CONVERT TO HIRES JPEG WITH: https://pdftojpg.me/

fldr = 'irri200' #sys.argv[1]
chrom = '01' #str(sys.argv[2]).zfill(2)
targets = ['lo', 'hi']
data = pd.read_csv(f"./irri200_GL_data.csv", header=0)

labs = {'subcont':'Sub-continent', 'africa':'Africa', 'china':'China', 'indo':'Sundaland', 'se.asia':'SE Asia', 'lat.am':'Americas', 'europe':'Europe', 'australia':'Australia'}
colDik = {'subcont':'crimson', 'africa':'green', 'china':'orange', 'indo':'dodgerblue', 'se.asia':'purple', 'lat.am':'black', 'europe':'gray', 'australia':'yellow'}
colDik2 = {'gray': (140, 140, 140), 'yellow': (0, 255, 255), 'purple': (128, 0, 128), 'orange': (0, 128, 255),'dkblue': (200, 0, 0),
           'dodgerblue': (255, 100, 100), 'black': (1, 0, 0), 'green': (0, 127, 0), 'crimson': (60, 20, 220), 'dkred': (0, 0, 150)}

types = ['Landrace', 'RD variety']
datGrps = data.groupby(['REGION','VAR']).size().reset_index().rename(columns={0:'cnt'})
datGrps['id'] = datGrps['REGION'] + '_' + datGrps['VAR'] #pivot Tab

#fnc
def writeText(im, blc, s, col, fnt, phrase, l):
    '''write text on dummy array same size as im'''
    font = fnt
    fontScale = s
    lineType = l #int(im.shape[0] / 250)

    txt = np.zeros(im.shape, np.uint8) # Create a black image - NB pointerlen determiness dist from pointer!!
    bottomLeftCornerOfText = blc
    fontColor = colDik2.get(col)
    cv2.putText(txt, phrase, bottomLeftCornerOfText, font, fontScale, fontColor, lineType)

    txtposns = np.where(txt != 0) #get txt posns in txt (i.e. that aren't white) and overwrite only these on im 
    im[txtposns[0], txtposns[1], :] = fontColor

    return im

for target in targets:
    df = pd.DataFrame()
    files = glob.glob1(fldr,f"*{chrom}g*{target}_assocSNPs.csv") #donors OR assocSNPs
    for file in files:
        potes = pd.read_csv(f"{fldr}/{file}", header=0).drop_duplicates()
        potes = potes[potes.type.isin(types)] #uncomment if using *assocSNPs.csv
        if potes.empty: continue
        #haps = pd.read_csv(f"{fldr}/{file.replace(f'_potential_{target}_donors', '.haps')}", header=0) #USE IF donors
        haps = pd.read_csv(f"{fldr}/{file.replace(f'_{target}_assocSNPs', '.haps')}", header=0) #USE IF assocSNPs
        haps['hap'] = haps.iloc[:, 1:-1].agg(''.join, axis=1)
        for al in potes.alleles.unique():
            haplo = ''.join(al.split('; '))
            samps = haps.acc[haps.hap == haplo]
            df = pd.concat([df, data[data.acc.isin(samps)]])

    #mk pivot & mge w prev pivot
    groups = df.groupby(['REGION','VAR']).size().reset_index().rename(columns={0:'cnt'})
    groups['id'] = groups['REGION'] + '_' + groups['VAR']
    groups = pd.merge(groups, datGrps, on='id')
    groups['adjusted'] = groups['cnt_x']/groups['cnt_y'] #adjust counts to reflect haps per accession

    #draw pies - segments are proportional to nReads per region
    for var in groups.VAR_x.unique():
        df2 = groups[groups.VAR_x == var]
        sizes = df2.cnt_x.tolist() #orig: df2.adjusted.tolist()
        cols, labels = [], [] #labels are nReads per accession
        [cols.append(colDik.get(x)) for x in df2.REGION_x.tolist()]
        [labels.append(str(round(df2.adjusted[df2.REGION_x == qw].iloc[0], 1))) for qw in df2.REGION_x.tolist()]

        fig1, ax1 = plt.subplots()
        ax1.pie(sizes, startangle=90, colors=cols, labels = labels)
        ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
        fig1.savefig(f"{var}_{target}_pie.png", dpi = 600)    

    #make panel plots
    pies = glob.glob1(".", f"*{target}_pie.png")
    crops = [] #crops = [variety type, cropped pie]

    for pie in pies:
        pic = cv2.imread(pie)[100:2600, 700:3200, :]
        crops.append([pie.split('_')[0], cv2.resize(pic, (0,0), fx=0.8, fy=0.8)]) 

    im = np.full((2400 * (1 + (int(len(crops)/2))), 3600, 3), 255, dtype=np.int16) #full white rectangle

    #Make legend
    Y = 1200 - ((75 * groups.REGION_x.unique().shape[0]) - 150) #should be good upto ca. 16 regions
    for lab in labs:
        im = writeText(im, (200, Y), 3, 'black', cv2.FONT_HERSHEY_SIMPLEX, labs.get(lab), 10)
        im[Y-100:Y+30, 1000:1130, :] = colDik2.get(colDik.get(lab))
        Y += 150

    #make legend border (region names max 15 chrs)
    Y2 = 1050 - ((75 * groups.REGION_x.unique().shape[0]) - 150)
    im[Y2:Y2+30, 100:1250, :] = colDik2.get('black')
    im[Y-100:Y-70, 100:1250, :] = colDik2.get('black')
    im[Y2:Y-70, 100:130, :] = colDik2.get('black')
    im[Y2:Y-70, 1220:1250, :] = colDik2.get('black')

    #pie panel plots
    y, x, z = 1200, 2700, -1 #pie position plus a functional var
    for crop in crops: #add pies
        txt = f"{crop[0]}: {groups.cnt_x[groups.VAR_x == crop[0]].sum()} accs."
        im = writeText(im, (x-700, y-1000), 4, 'black', cv2.FONT_HERSHEY_SIMPLEX, txt, 15)
        yDispl, xDispl = int(crop[1].shape[0] / 2), int(crop[1].shape[1] / 2)
        y1 = y - yDispl
        y2 = y1 + crop[1].shape[0]
        x1 = x - xDispl
        x2 = x1 + crop[1].shape[1]

        im[y1:y2, x1:x2, :] = crop[1]
        if z == 1: x += 1800
        else:
            y += 2400
            x = 900 #back to lhs x-position
        z *= -1

    #output figures for low & hi data
    #NB if > 7200 pixels long, will divide into seperate figures
    size, n = im.shape[0], 1
    start, stop, incr = 0, 7200, 7200
    while size + incr > stop:
        cv2.imwrite(f"Chrom{chrom}_PieChart_{target}_image{n}.png", im[start:min(stop, size), :, :])
        start += incr
        stop += incr
        n += 1

    for png in glob.glob1(".", f"*{target}_pie.png"): os.remove(png)

