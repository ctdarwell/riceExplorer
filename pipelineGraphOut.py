# sys.argv[1-3]: work dir; gene annotations; genes of interest list
# OPTIONAL COHORT sys.argv[4-6]: chrom; prefix; known gene position in CSV (integer)
import pandas as pd, glob, numpy as np, cv2, sys
import itertools
from PIL import ImageFont, ImageDraw, Image, ImageOps

#main params
wdir = sys.argv[1] #eg 'fullChromAnal'
MSU = sys.argv[2] # eg 'msu_ver7_simple.csv' 
GENE_LIST = sys.argv[3] # #eg msu_GL.csv
sign = sys.argv[4] #'LOC_Os' # 

fig1 = sys.argv[5].replace('_', ' ').split(',') #NFV types  sys.argv[5]
fig2 = sys.argv[6].replace('_', ' ').split(',') #EC types sys.argv[6]

#other params
s, loc_wd, pntr_wd, pntr_len, off, rh_margin = 10, 15, 6, 150, 10, 2 #size of chrom, graphical locus marker width, pointer width, pointer length, txt offset, indent size
FILES = glob.glob1(wdir, f'{sign}*_assocSNPs.csv')
colDik = {'white': (255, 255, 255), 'orng': (255, 140, 0), 'dkred': (0, 0, 150), 'blue': (0, 150, 255), 'red': (0, 0, 255), 'black': (1, 0, 0), 'green': (0, 200, 0)}
try: prfx = f" - {sys.argv[14]}_"
except: prfx = ""

#read files
msu = pd.read_csv(MSU, header=None)

col_dik = {int(sys.argv[7]) - 1: 'chr', int(sys.argv[8]) - 1: 'locus', int(sys.argv[9]) - 1: 'start', int(sys.argv[10]) - 1: 'stop'}
columns = [col_dik.get(i) if i in col_dik else f"Col{i+1}" for i in range(msu.shape[1])]
msu.columns = columns

try: gene_list = pd.read_csv(f"{GENE_LIST}", header=None)
except: gene_list = None

def main():
    '''Main pgm control function'''
    chroms = [x.split(sign)[1][:len(sys.argv[11])] for x in FILES] #

    try: chroms = [sys.argv[15].zfill(2)] #single chrom for custom figure
    except: pass

    for chrom in sorted(set(chroms)):
        #mk blank img + chrom
        im = np.full((3600, 3600, 3), 255, dtype=np.int16) #full white rectangle #ALT!!!!!!!!!!!!!
        h, l = [int((im.shape[0]/3) - 50), int((im.shape[0]/3) + 50)], [int(im.shape[1]/s), int(im.shape[1] - (im.shape[1]/(s / rh_margin)))]
        im[h[0]:h[1], l[0]:l[1], :2] = 0 #colour bar
        for p in range(1, 26): im[h[0] + p:h[1] - p, l[0] - p:l[1] + p, :2] = 0 #colour bar ends
        h2 = [int(h[0] + (1550)), int(h[1] + (1550))]
        im[h2[0]:h2[1], l[0]:l[1], :] = (0,255,0) #colour bar
        for p in range(1, 26): im[h2[0] + p:h2[1] - p, l[0] - p:l[1] + p, :] = (0,255,0) #colour bar ends

        #get data from output files
        genes, snp_midpnts, catalog = getData(chrom)
        
        #establish chunk of chrom to annotate (eg if snps range: 6,764,235 to 38,953,002 will illustrate 5-40Mb)
        start_end = [int(f"{str(snp_midpnts[0]-1000000)[:-6]}000000"), int(f"{str(snp_midpnts[-1])[:-6]}000000") + 2000000]
        
        #mark known important genes
        try: im = mkGene(im, start_end, h, h2, l, chrom)
        except: pass

        #accession types of interest only        
        cat1 = catalog[catalog.type.isin(fig1)]
        cat2 = catalog[catalog.type.isin(fig2)]
        intersect = cat2.locus.unique()[np.in1d(cat2.locus.unique(), cat1.locus.unique())]

        #convert SNPs to pixels in range of array - and mark genes on chromosome bar
        pixels1, genes1 = getPixels(im, genes, snp_midpnts, start_end, l, h, loc_wd, cat1)
        pixels2, genes2 = getPixels(im, genes, snp_midpnts, start_end, l, h2, loc_wd, cat2)

        #separate hi vs. lo value genotypes
        hipix1 = pixels1[np.where(cat1.locus[cat1.hilo == 'hi'].unique()[:, None] == genes1.ravel())[1]]
        hignz1 = genes1[np.where(cat1.locus[cat1.hilo == 'hi'].unique()[:, None] == genes1.ravel())[1]]
        lopix1 = pixels1[np.where(cat1.locus[cat1.hilo == 'lo'].unique()[:, None] == genes1.ravel())[1]]        
        lognz1 = genes1[np.where(cat1.locus[cat1.hilo == 'lo'].unique()[:, None] == genes1.ravel())[1]]        
        
        hipix2 = pixels2[np.where(cat2.locus[cat2.hilo == 'hi'].unique()[:, None] == genes2.ravel())[1]]
        hignz2 = genes2[np.where(cat2.locus[cat2.hilo == 'hi'].unique()[:, None] == genes2.ravel())[1]]
        lopix2 = pixels2[np.where(cat2.locus[cat2.hilo == 'lo'].unique()[:, None] == genes2.ravel())[1]]
        lognz2 = genes2[np.where(cat2.locus[cat2.hilo == 'lo'].unique()[:, None] == genes2.ravel())[1]]
        
        #unsquish any overlapping labels - COULD JUST SEND im.shape HERE!!!!!
        hiEnds1 = unSquishTxt(im, hipix1) #labs1 & labs2 denote whether squished labels could not be separated....
        loEnds1 = unSquishTxt(im, lopix1) #...will result in no labels being added
        hiEnds2 = unSquishTxt(im, hipix2) #labs1 & labs2 denote whether squished labels could not be separated....
        loEnds2 = unSquishTxt(im, lopix2)

        #draw pointers
        im = pntrs(im, hipix1, hiEnds1, pntr_wd, h[0] - 50, pntr_len * -1, -1, 1)
        im = pntrs(im, lopix1, loEnds1, pntr_wd, h[1] + 50, pntr_len, 1, -1)
        im = pntrs(im, hipix2, hiEnds2, pntr_wd, h2[0] - 50, pntr_len * -1, -1, 1)
        im = pntrs(im, lopix2, loEnds2, pntr_wd, h2[1] + 50, pntr_len, 1, -1)

        PIL_image = Image.fromarray(np.uint8(im)).convert('RGB')
        
        #TO BE TESTED inc OUTPUT FOLDER!!!!!!!!!!!!!!!!!!!!!!!!!!
        outputNames, c = ['hignz1', 'hignz2', 'lognz1', 'lognz2'], 0
        for output in [hignz1, hignz2, lognz1, lognz2]:
            outputDF = pd.DataFrame(output)
            outputDF.columns = ['loci']
            outputDF.to_csv(f"{wdir}/LOC_Os{chrom}_{outputNames[c]}.csv", index=False)
            c += 1

        #write gene names (!!genes labels may be omitted if labels are squashed)
        if hignz1.shape[0] < 63: PIL_image = rotateText(hiEnds1, hignz1, 90, 30, 575, intersect, PIL_image)
        if lognz1.shape[0] < 63: PIL_image = rotateText(loEnds1, lognz1, 270, 30, 1500, intersect, PIL_image)
        if hignz2.shape[0] < 63: PIL_image = rotateText(hiEnds2, hignz2, 90, 30, 2125, intersect, PIL_image)
        if lognz2.shape[0] < 63: PIL_image = rotateText(loEnds2, lognz2, 270, 30, 3050, intersect, PIL_image)

        draw = ImageDraw.Draw(PIL_image)

        font = ImageFont.truetype(sys.argv[12], 90) # use a truetype font
        draw.text((700, 150), f"CHROMOSOME {int(chrom)}: {int(start_end[0]/1000000)}-{int(start_end[1]/1000000)}Mb{prfx[:-1]}", font=font, fill=(0, 0, 150, 0)) #I think the last value in fill is opacity (after rgb)
        font = ImageFont.truetype(sys.argv[13], 60)
        dik = {(3000, 850): 'NFVs: High-value', (3000, 925): 'genotypes', (3000, 1400): 'NFVs: Low-value', (3000, 1475): 'genotypes', (3000, 2400): 'EVs: High-value', (3000, 2475): 'genotypes', (3000, 2950): 'EVs: Low-value', (3000, 3025): 'genotypes'}
        for d in dik: draw.text(d, dik.get(d), font=font, fill=(0, 0, 0, 0))

        PIL_image.save(f"{wdir}/{prfx[3:]}Chrom{chrom}_potential_donors.jpg", dpi=(600, 600))
        pdf = PIL_image.convert('RGB')
        pdf.save(f"{wdir}/{prfx[3:]}Chrom{chrom}_potential_donors.pdf", dpi=(600, 600))

#MODIFIED TO MARK ANY KNOWN GENES WITH ANY FROM ANAL [COULD MOFIFY COLOUR ON OUTPUT]
#NB THIS FUNC ONLY CALLED IF ADDITIONAL SYSARGS ENTERED
def mkGene(im, start_end, h, h2, l, chrom):
    tmp = [(255, 140, 0), (0, 150, 255), (0, 255, 255), (135, 206, 235), (0, 0, 128), (0, 139, 139), (70, 130, 180), (0, 90, 156), (0, 128, 255), (0, 127, 0), (255, 100, 100), (140, 140, 140)]
    cnt, step = 0, 0
    for gene in gene_list.index:
        try:
            if gene == int(sys.argv[16]) - 1:
                known = int(sys.argv[16])
                step, cnt = 1, 1
                for indx in range(known - 1, gene_list.shape[0]):tmp.append((gene_list.iloc[indx,1], gene_list.iloc[indx,2], gene_list.iloc[indx,3])) 
        except: pass
        if f"{sign}{chrom}" not in gene_list.iloc[gene, 0]: continue
        mu = int(np.mean([float(msu.stop[msu.locus == gene_list.iloc[gene, 0]].iloc[0]), float(msu.start[msu.locus == gene_list.iloc[gene, 0]].iloc[0])]))
        pix = int(l[0] + ((mu - start_end[0]) / (start_end[1] - start_end[0]) * (l[1] - l[0]))) #calc relative posn on bar
        if pix > l[1] or pix < l[0]: continue #if not in range
        im[h[1]:h[1] + 25, pix - loc_wd:pix + loc_wd, :] = tmp[cnt] # :2] =0 OR :] = tmp[cnt]
        im[h[0] - 25:h[0], pix - loc_wd:pix + loc_wd, :] = tmp[cnt] # :2] =0 OR :] = tmp[cnt]
        im[h2[1]:h2[1] + 25, pix - loc_wd:pix + loc_wd, :] = tmp[cnt] # :2] =0 OR :] = tmp[cnt]
        im[h2[0] - 25:h2[0], pix - loc_wd:pix + loc_wd, :] = tmp[cnt] # :2] =0 OR :] = tmp[cnt] 
        cnt+=step

    return im

def getPixels(im, genes, snp_midpnts, start_end, l, h, loc_wd, cat):
    indxs = np.where(cat.locus.unique()[:, None] == genes.ravel())[1]

    dik = {'lo': [int(np.mean(h)), h[1]], 'hi': [h[0], int(np.mean(h))]}
    dik2 = {'lo': [int(np.mean(h)), int(np.mean(h)) + 25], 'hi': [int(np.mean(h)) - 25, int(np.mean(h))]}

    pixels, loci = [], []
    for indx in indxs:
        pix = int(l[0] + ((snp_midpnts[indx] - start_end[0]) / (start_end[1] - start_end[0]) * (l[1] - l[0]))) #calc relative posn on bar
        for hilo in cat.hilo[cat.locus == genes[indx]].unique():
            im[dik.get(hilo)[0]:dik.get(hilo)[1], pix - loc_wd:pix + loc_wd, :] = 255 #mark genes with black dash
            im[dik2.get(hilo)[0]:dik2.get(hilo)[1], pix - loc_wd:pix + loc_wd, :] = 0
        if pix in pixels: continue
        pixels.append(pix), loci.append(genes[indx])

    return np.array(pixels), np.array(loci)

def getData(chrom):
    snp_midpnts, genes, catalog = [], [], pd.DataFrame() #SNPs RM!
    files = sorted(glob.glob1(wdir, f'{sign}{chrom}*_assocSNPs.csv')) #both hi & lo val files
    for f in files:
        file = pd.read_csv(f"{wdir}/{f}", header=0)
        df = file[['acc','type','locus']]
        df['hilo'] = f.split('_')[-2] #add column with hi or low designation
        df = df[df.type.isin(fig1 + fig2)] #ignore file if acc type not of interest
        if df.empty: continue
        catalog = pd.concat([catalog, df])
        genes.append(df.locus.unique()[0])
        snp_midpnts.append(int(np.mean([float(msu.stop[msu.locus == df.locus.unique()[0]].iloc[0]), float(msu.start[msu.locus == df.locus.unique()[0]].iloc[0])])))

    return np.array(genes), np.array(snp_midpnts), catalog.drop_duplicates().reset_index(drop=True)

#annotate markers
def rotateText(pointerLocs, genes, angle, off, d, intersect, PIL_image):
    '''write text on dummy array same size as im
    then rotate before applying to im with appropriate colors'''

    xaposns = dict(zip(genes, pointerLocs))
    font = ImageFont.truetype(sys.argv[13], 40) # use a truetype font

    for gene in genes:
        if gene in intersect: fontColor = 255
        else: fontColor = 0
        txt = Image.new('L', (340, 50), color=0) #mk text box
        imd = ImageDraw.Draw(txt)
        imd.text((0, 0), gene,  font=font, fill=255) #0=black
        w = txt.rotate(angle,  expand=1)
        
        PIL_image.paste(ImageOps.colorize(w, (0, 128, 0), (0, 0, fontColor)), (xaposns.get(gene) - off, d),  w) #w, ?, colour, posn

    return PIL_image

def writeText(im, blc, s, col, fnt, phrase):
    '''write text on dummy array same size as im'''
    font = fnt
    fontScale = s
    lineType = int(im.shape[0] / 250)

    txt = np.zeros(im.shape, np.uint8) # Create a black image - NB pointerlen determiness dist from pointer!!
    bottomLeftCornerOfText = blc
    fontColor = colDik.get(col)
    cv2.putText(txt, phrase, bottomLeftCornerOfText, 
                font, fontScale, fontColor, lineType)

    txtposns = np.where(txt != 0) #get txt posns in txt (i.e. that aren't white) and overwrite only these on im 
    im[txtposns[0], txtposns[1], :] = fontColor

    return im

def unSquishTxt(im, pointerLocs):
    '''check there is enough space between the labels
    re-position if too squished'''
    gap = 45 #int(im.shape[0] / 30)
    dists, probs, cnt = [], [[]], 0
    [dists.append(x) for x in pointerLocs] #pointerLocs is XYs
    dists.insert(0, 100), dists.append(im.shape[1] - 100) #ADDED EXPT LINE - imaginery locations near the limits of the im so calcs below have values to check if they're near something (i.e. the list indexes don't get out of bounds) - a bit fudgey
    for dist in dists[:-1]:
        #make list of problematic posns
        if dists[dists.index(dist) + 1] - dist < gap: probs[cnt].extend([dist, dists[dists.index(dist) + 1]])
        else:
            probs.append([])
            cnt += 1

    #house keeping
    while [] in probs: probs.remove([])
    if probs == []: return np.array(pointerLocs)
    
    swch, cntr = [1], 0
    while swch != []:
        cntr += 1
        if cntr > 100: return np.linspace(50, 2880, len(pointerLocs)).astype(int) #np.array(pointerLocs)
        swch = []
        #merge any clusters of problems together if they contain a common problem posn
        for it in itertools.combinations(probs, 2):
            if list(set(it[0]) & set(it[1])) != []:
                try: probs.remove(probs[probs.index(it[1])])
                except: return np.linspace(50, 2880, len(pointerLocs)).astype(int)
                try: probs.remove(probs[probs.index(it[0])])
                except: return np.linspace(50, 2880, len(pointerLocs)).astype(int)
                probs.insert(0, it[0] + it[1])
                
        for prob in probs:
            #calculate new gaps btwn txt (i.e. x coords)
            setprob = sorted(list(set(prob)))
            mu = int(np.mean([setprob[0], setprob[-1]]))
            init = mu - ((len(setprob) - 1) * int(gap / 2))
            for i in setprob[:-1]: init += gap
            #check new posns aren't problematic - otherwise back to top of loop and redo w new posns
            if mu - ((len(setprob) - 1) * int(gap / 2)) - dists[dists.index(setprob[0]) - 1] < gap:
                probs[probs.index(prob)].append(dists[dists.index(setprob[0]) - 1])
                swch.append(1) #any issues mean swch != []: thus, goes back to top of while loop
            
            try:
                if dists[dists.index(setprob[-1]) + 1] - init < gap:
                    probs[probs.index(prob)].append(dists[dists.index(setprob[-1]) + 1])
                    swch.append(1)
            except: return np.linspace(50, 2880, len(pointerLocs)).astype(int)

    #assign new gaps
    new_arr = np.array(pointerLocs).copy()
    for prob in probs:
        setprob = sorted(list(set(prob))) #recalc these
        mu = int(np.mean([setprob[0], setprob[-1]]))
        init = mu - ((len(setprob) - 1) * int(gap / 2))
        for i in setprob: #change vals in pointerLocs
            new_arr[np.where(new_arr == i)] = init
            init += gap

    #double check no overlapping of borders - return default if so
    if new_arr[-1] > 2880: new_arr = new_arr - (new_arr[-1] - 2880)
    if new_arr[0] < 50: new_arr = new_arr + (50 - new_arr[0])
    if new_arr[-1] > 2880: return np.linspace(50, 2880, len(pointerLocs)).astype(int)

    return new_arr

#draw pointers
def pntrs(im, hilo, ends, pntr_wd, start, y, z, q):
    for loc in hilo:
        end = ((loc - ends[np.where(hilo == loc)][0]) / y) * q
        for i in range(0, y, z):
            im[start + i, int(loc)-pntr_wd:int(loc)+pntr_wd] = colDik.get('black')
            loc += end

    return im

if __name__ == '__main__': main()

