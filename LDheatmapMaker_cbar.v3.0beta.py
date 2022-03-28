#input file generated in vcftools
#1. automate LOC_Os upload - see pipelineGraphOut_hiRes.faoGithub.py
#2. improve text
from PIL import ImageFont, ImageDraw, Image, ImageOps
import pandas as pd, numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import cv2, itertools, glob, sys
from scipy import ndimage


file = sys.argv[1] #'*LD.csv'
workfldr = sys.argv[2]
chrom = file[3:5]
HAPS =  sys.argv[3] #'LOC_Os07_hignz1.csv'

#load full gene list and assigne columns
msu = pd.read_csv(sys.argv[4], header=0) #INPUT
msu.columns = [f'col{x+1}' for x in range(msu.shape[1])]
msu_cols = msu.columns.tolist()
msu_cols[int(sys.argv[7]) - 1] = 'locus'
msu_cols[int(sys.argv[8]) - 1] = 'start'
msu_cols[int(sys.argv[9]) - 1] = 'stop'
msu.columns = msu_cols

genes = sorted(pd.read_csv(f"{workfldr}/{HAPS}", header=0).loci.tolist())
search = sys.argv[5] #'LOC_Os' #INPUT

start = msu.start[msu.locus.isin(genes)].min() - 500000
end = msu.stop[msu.locus.isin(genes)].max() + 500000
region = [start, end] # end + 2000000 may work! [start, end] OR None; chr region shown on img

try:
    known = pd.read_csv(sys.argv[6], header=None) #TRY/EXCEPT
    interest = known[0][known[0].str.startswith(f"{search}{chrom}")].tolist()
    toadd = sorted(msu.locus[(msu.start > start) & (msu.stop < end) & (msu.locus.isin(interest))].tolist())
    genes.extend(toadd)
except: pass

genes = sorted(set(genes)) #NB id'd genes will be overwritten by any gene of interest
colours = ['dkred'] * len(genes)

try:
    for _ in toadd: colours[genes.index(_)] = 'black'
except: pass

markers_to_highlight = []
for gene in genes:
    markers_to_highlight.append(int((msu.stop[msu.locus == gene] - msu.start[msu.locus == gene]) / 2) + int(msu.start[msu.locus == gene]))

#clusters (regions of interest) of 3 or more markers at 200000bp
rois = []
grps = np.split(markers_to_highlight, np.where(np.diff(markers_to_highlight) > 200000)[0] + 1)
for grp in grps:
    if grp.size > 2: rois.append([max(region[0], grp[0] - 95000), min(grp[-1] + 95000, region[1])])

colDik = {'white': (255, 255, 255), 'blue': (200, 0, 0), 'red': (0, 0, 255), 'black': (1, 0, 0), 'green': (0, 200, 0), 'dkred': (0, 0, 150)}
roicol = colDik.get('dkred') #can change region of interest colour
pointercols = [colDik.get(c) for c in colours]
pointerdik = dict(zip(markers_to_highlight, pointercols))

def main():
    '''
    Main pgm control function
    '''
    if (len(genes) != len(markers_to_highlight)) or (len(genes) != len(pointercols)) or (len(markers_to_highlight) != len(pointercols)):
        print('RUN STOPPED: variables "genes", "markers_to_highlight", and "pointercols" are not all of equal length')
        return
    pointers = markers_to_highlight.copy()
    dat = pd.read_csv(file, header=0)
    print(f"Read: {file}; size = {dat.shape}")
    im, cols = mkHMP(dat)
    starty, startx, endx = getXYcoords(im)
    snpmap = [[startx, endx], region]
    maplen = snpmap[1][1] - snpmap[1][0]
    loci = dat.columns[cols[:-1]]
    rois, pointers = checkMarkers(loci, pointers)
    pointerlen, pointerWid, pointerLocs = calcPointers(im, pointerdik, snpmap, maplen, loci, startx, endx, starty, pointers)
    im = mkROIs(im, starty, loci, startx, endx, rois)
    pointerLocs = unSquishTxt(im, pointerLocs, snpmap) #MIGHT GIVE PROBS AT SOME POINT!!

        
    im = mkPointers(im, pointerdik, loci, startx, endx, starty, pointerLocs, pointers)

    PIL_image = Image.fromarray(cv2.cvtColor(im, cv2.COLOR_BGR2RGB))
    PIL_image = rotateText(np.array(pointerLocs), genes, 90, 30, 1080, colours, PIL_image)

    #CAN USE: print(f"\nFiles = {HAPS.replace(workfldr, f'{workfldr}/outFiles').replace('.csv', '')}_{int(region[0]/1000000)}-{int(region[1]/1000000)}Mb_LD.*")
    print(f"\nFiles = {HAPS.replace('.csv', '')}_{int(region[0]/1000000)}-{int(region[1]/1000000)}Mb_LD.*")
    PIL_image.save(f"{HAPS.replace('.csv', '')}_{int(region[0]/1000000)}-{int(region[1]/1000000)}Mb_LD.png", dpi=(600, 600))
    PIL_image.save(f"{HAPS.replace('.csv', '')}_{int(region[0]/1000000)}-{int(region[1]/1000000)}Mb_LD.pdf", dpi=(600, 600))

def unSquishTxt(im, pointerLocs, snpmap):
    '''check there is enough space between the labels
    re-position if too squished'''
    gap = int(im.shape[0] / 60) #orig val = 50
    probs, cnt = [[]], 0
    dists = [x[0] for x in pointerLocs]
    dists.insert(0, 100), dists.append(im.shape[1] - 100) #ADDED EXPT LINE - imaginery locations near the limits of the im so calcs below have values to check if they're near something (i.e. the list indexes don't get out of bounds) - a bit fudgey
    for dist in dists[:-1]:
        #make list of problematic posns
        if dists[dists.index(dist) + 1] - dist < gap:
            probs[cnt].extend([dist, dists[dists.index(dist) + 1]])
        else:
            probs.append([])
            cnt += 1

    #house keeping
    while [] in probs: probs.remove([])
    if probs == []: return pointerLocs
    
    swch, cntr = [1], 0
    while swch != []:
        cntr += 1
        if cntr > 3: #give up after 3 goes and give default spread of pointer positions
            start = max(100, int(im.shape[0] / 2) - (55 * len(pointerLocs) / 2))
            end = min(int(im.shape[0] / 2) + (55 * len(pointerLocs) / 2), im.shape[0] - 100)
            new_coords = np.linspace(start, end, len(pointerLocs)).astype(int)
            return list(zip(new_coords, np.array(pointerLocs)[:,1]))

        swch = []
        #merge any clusters of problems together if they contain a common problem posn
        for it in itertools.combinations(probs, 2):
            if list(set(it[0]) & set(it[1])) != []:
                try: probs.remove(probs[probs.index(it[1])])
                except: pass
                try: probs.remove(probs[probs.index(it[0])])
                except: pass
                probs.insert(0, it[0] + it[1])
                
        for prob in probs:
            #calculate new gaps btwn txt (i.e. x coords)
            setprob = sorted(list(set(prob)))
            mu = int(np.mean([setprob[0], setprob[-1]])) #prob midpoint
            init = mu - ((len(setprob) - 1) * int(gap / 2)) #new start point?
            for i in setprob[:-1]: init += gap #new end point?
            #check new posns aren't problematic - otherwise back to top of loop and redo w new posns
            if mu - ((len(setprob) - 1) * int(gap / 2)) - dists[dists.index(setprob[0]) - 1] < gap:
                probs[probs.index(prob)].append(dists[dists.index(setprob[0]) - 1])
                swch.append(1) #any issues mean swch != []: thus, goes back to top of while loop
            if dists[dists.index(setprob[-1]) + 1] - init < gap:
                probs[probs.index(prob)].append(dists[dists.index(setprob[-1]) + 1])
                swch.append(1)
                #break #- may be helpful for trouble shooting

    #assign new gaps
    for prob in probs:
        setprob = sorted(list(set(prob))) #recalc these
        mu = int(np.mean([setprob[0], setprob[-1]]))
        init = mu - ((len(setprob) - 1) * int(gap / 2))
        for i in setprob: #change vals in pointerLocs
            pointerLocs[pointerLocs.index([i, pointerLocs[0][1]])][0] = init
            init += gap

    return pointerLocs  

#annotate markers - orig from pipelineGraphOut_hiRes.faoGithub.py
def rotateText(pointerLocs, genes, angle, off, d, colours, PIL_image):
    '''write text on dummy array same size as im
    then rotate before applying to im with appropriate colors'''
    xaposns = dict(zip(genes, pointerLocs[:, 0] + 8))
    font = ImageFont.truetype(sys.argv[10], 40) # use a truetype font
    
    for gene in genes:
        fontColor = colDik.get(colours[genes.index(gene)])[::-1] #colours reversed for pil image
        txt = Image.new('L', (340, 50), color=0) #mk text box
        imd = ImageDraw.Draw(txt)
        imd.text((0, 0), gene,  font=font, fill=255) #0=black
        w = txt.rotate(angle,  expand=1)
        
        PIL_image.paste(ImageOps.colorize(w, (0, 128, 0), fontColor), (xaposns.get(gene) - off, d),  w) #w, ?, colour, posn

    return PIL_image

def mkPointers(im, pointerdik, loci, startx, endx, starty, pointerLocs, pointers):
    '''calculate start and end point of pointers and draw.
    Bottom of pointer is at a SNP featured in the data (maybe recalculated in 'checkMarkers' if non-present marker given).
    Top of pointer is original given position, irrespective of whether the position has a SNP in the data or not'''    
    pointerlen, pointerGap, pointerWid = int(im.shape[0] / 20), int(im.shape[0] / 280), int(im.shape[0] / 350)
    for p in pointers:
        try: p0, p1 = p[0], p[1] #if pointer has recalculated value alongside original that has no representative marker
        except: p0, p1 = p, p
        pointerTOP = pointerLocs[pointers.index(p)][0] # int(startx + (pointprop * (endx - startx))) #start point is proportional to end/start of map
        pointerBOT = startx + int(np.where(loci == str(p1))[0][0] / loci.shape[0] * (endx - startx)) #bottom is proportional according to no. of markers
        incr = (pointerBOT - pointerTOP) / pointerlen #calc x increment per y pixel
        for ypos in range(starty - pointerlen - pointerGap, starty - pointerGap): #draw the lines
            im[ypos:ypos + 1, int(pointerTOP + incr): int(pointerTOP + incr + pointerWid)] = pointerdik.get(p0)
            pointerTOP += incr

    return im

def mkHMP(dat):
    '''make the heat map'''
    print('\nHeatmap to be annotated [MAY BE BLANK!!]')
    cols = []
    if region:
        for col in dat.columns[1:]: #id df columns in chosen region (also useful later on)
            if float(col) > region[0]: cols.append(dat.columns.tolist().index(col))
            if float(col) > region[1]: break
    else: cols = [0, len(dat.columns) - 1]
    arr = np.array(dat.iloc[cols[0]:cols[-1], cols[0] + 2:cols[-1] + 2]) #extract stipulated region
    mask = np.zeros_like(arr)
    mask[np.tril_indices_from(mask)] = True #something to do with stating upper tri as lower
    
    #create a heatmap - cbar=True gives colorbar (fao extracting colour bar - heatmap is too small)
    fig = plt.figure(1)
    with sns.axes_style("white"):
        sns.heatmap(arr, xticklabels=False, cbar=True, 
                    cbar_kws={"orientation": "horizontal", "shrink": 0.18},
                    yticklabels=False, mask=mask, square=True,  cmap="YlOrRd")

    fig.savefig("tmp.png", dpi = 450) #save so can open with cv2
    plt.close(fig)

    #extract color bar
    im = cv2.imread('tmp.png')
    bar, ys, xs = getBarCoords(im)
    
    #create 2nd heatmap - cbar=False gives no colorbar (the heatmap is now bigger)
    fig = plt.figure(1)
    with sns.axes_style("white"):
        sns.heatmap(arr, xticklabels=False, cbar=False, 
                    yticklabels=False, mask=mask, square=True,  cmap="YlOrRd")

    fig.savefig("tmp.png", dpi = 450)
    plt.close(fig)

    #rotate heatmap and add colour bar
    im = ndimage.rotate(cv2.imread("tmp.png"), 225, mode='nearest') #rotate by 225 deg
    y, x = int(im.shape[0] * .75), int(im.shape[0] * .565) #calc where to put bar
    im[y:y+bar.shape[0], x:x+bar.shape[1], :] = bar #overlay bar

    return im, cols

def getBarCoords(im2): #this could feasibly be mgd w getXYcoords as a single func
    '''similar to "getXYcoords" but finds outer edges of the colour bar
    then returns its coords'''
    l1 = []
    for starty2 in range(im2.shape[0] - 1, 0, -1): #pixel by pixel up image's centre til non-white pixel found 
        if im2[starty2, im2.shape[1]//2, 0] != 255 or im2[starty2, im2.shape[1]//2, 1] != 255 or im2[starty2, im2.shape[1]//2, 2] != 255:
            l1.append(starty2)

    for l in l1: #evaluate top & bottom border of bar - list of non-white pixels need splitting into bar vs. heatmap
        if l - l1[l1.index(l) + 1] > im2.shape[0] / 8: break #finds where the bar ends (going upwards) - i.e. the gap btwn the bar and the actual heat map (there are gaps of white within the bar so this is a bit of fiddling to find the real gap)
    ys = [l1[0] + int(im2.shape[0] / 80), l - int(im2.shape[0] / 80)] #top & bottom of bar

    for startx2 in range(im2.shape[1]): #find LH edge of bar using upper edge y-coord of bar
        if im2[l, startx2, 0] != 255 or im2[l, startx2, 1] != 255 or im2[l, startx2, 2] != 255:
            break

    xs = [startx2 - int(im2.shape[0] / 10), im2.shape[1] - startx2 + int(im2.shape[0] / 10)]
    bar = im2[ys[1]: ys[0], xs[0]:xs[1]] #coords of bar with a border

    return bar, ys, xs

def getXYcoords(im):
    '''identify start/end posns (i.e. margins) of heat map against blank background
    work down y-axis from midpoint of x (we know heat map is in centre of img
    from id'd y, work inwards to find LHS x value
    endx given by symmetry of shape'''
    for starty in range(im.shape[0]): #pixel by pixel down image's centre til non-white pixel found
        if im[starty, im.shape[1]//2, 0] != 255 or im[starty, im.shape[1]//2, 1] != 255 or im[starty, im.shape[1]//2, 2] != 255:
            break

    for startx in range(im.shape[1]): #find LH edge of heatmap using upper edge y-coord of bar
        if im[starty, startx, 0] != 255 or im[starty, startx, 1] != 255 or im[starty, startx, 2] != 255:
            break

    endx = im.shape[1] - startx #RH edge of heatmap

    return starty, startx, endx

def checkMarkers(loci, pointers):
    '''check stipulated markers are actually in the data -
    if not, for ROIs, take markers giving widest region within the stipulated boundaries
    for pointers, take nearest marker - but keep record of original input as this chsomal posn will be mapped'''

    newrois, newpointers = [], []
    mks = loci.astype(np.int64)
    for roi in rois:
        #if roi[0] not in mks: print(f"No recorded SNP at {roi[0]}!\nAltering region of interest values accordingly.")
        #if roi[1] not in mks: print(f"No recorded SNP at {roi[1]}!\nAltering region of interest values accordingly.")
        newrois.append([mks[np.where(mks >= roi[0])[0][0]], mks[np.where(mks > roi[1])[0][0] - 1]]) #get next inward SNPs from stipulated start-end of ROI
    print('\n')
    for pointer in pointers: 
        if pointer not in mks: #below will store two values of pointer as list within list of pointers (i.e. inputted value shows chromosomal posn while calculated nearest value used to mark on LD map) 
            #print(f"No recored SNP at {pointer}!\nAltering gene position to nearest marker")
            v1, v2 = abs(pointer - mks[np.where(mks > pointer)[0][0] - 1]), abs(pointer - mks[np.where(mks > pointer)[0][0]])
            if v1 < v2: newpointers.append([pointer, mks[np.where(mks > pointer)[0][0] - 1]])
            else: newpointers.append([pointer, mks[np.where(mks > pointer)[0][0]]])
        else: newpointers.append(pointer)

    return newrois, newpointers

def calcPointers(im, pointerdik, snpmap, maplen, loci, startx, endx, starty, pointers):
    '''calculate start and end point of pointers and draw.
    Bottom of pointer is at a SNP featured in the data (maybe recalculated in 'checkMarkers' if non-present marker given).
    Top of pointer is original given position, irrespective of whether the position has a SNP in the data or not'''
    pointerLocs = []
    pointerlen, pointerGap, pointerWid = int(im.shape[0] / 20), int(im.shape[0] / 280), int(im.shape[0] / 350)
    for p in pointers:
        try: p0 = p[0] #p0, p1 = p[0], p[1] #if pointer has recalculated value alongside original that has no representative marker
        except: continue #p0, p1 = p, p
        pointprop  = (p0 - snpmap[1][0]) / maplen 
        pointerTOP = int(startx + (pointprop * (endx - startx))) #start point is proportional to end/start of map
        pointerLocs.append([pointerTOP, starty - pointerlen - pointerGap]) #need a record as referenc for later drawing

    return pointerlen, pointerWid, pointerLocs

def mkROIs(im, starty, loci, startx, endx, rois):
    '''draw triangles marking regions of interest
    maybe recalculated in 'checkMarkers' if locations given don't have a SNP in the data'''
    for roi in rois: #again calc like pointerBOT
        roistart = startx + int(np.where(loci == str(roi[0]))[0][0] / loci.shape[0] * (endx - startx))
        roiend = startx + int(np.where(loci == str(roi[1]))[0][0] / loci.shape[0] * (endx - startx))
        lim1, lim2 = int(im.shape[0] / 350), int(im.shape[0] / 300) #these arbitrary numbers control thickness of triangles

        #draw a hollow triangle
        y2 = starty - 5 #top of triangle
        l, r, d = roistart, roiend, 1
        while d < lim1: #width of horizontal line of tri
            im[y2 + d -1: y2 + d, l:r] = roicol
            l += 1
            r -= 1
            d += 1
        while r - lim2 > l: #draw angled lines of tri with gaps between (ie don't overwrite heatmap)
            im[y2 + d -1: y2 + d, l:l + lim2] = roicol
            im[y2 + d -1: y2 + d, r - lim2:r] = roicol
            l += 1
            r -= 1
            d += 1
        while r > l: #make tip w/o gap btwn lines
            im[y2 + d -1: y2 + d, l:r] = roicol
            l += 1
            r -= 1
            d += 1

    return im


if __name__ == '__main__': main()
