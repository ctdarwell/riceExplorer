import glob, sys, numpy as np
import pandas as pd
from PIL import ImageFont, ImageDraw, Image, ImageOps

job = sys.argv[2]
length, width, w = 3000, 3600, 100 #w is dist of start of colored panel from image border

#get taxon names
qfile = f'sNMF_CL_v1.2/{sys.argv[1]}' #'wasps/run3.K8.Q'
vcf, fldr = f'{job}.allChroms.vcf', '.' # sys.argv[2], sys.argv[3] #'figwasp1617.vcf', '.'
for line in open(vcf):
    if line.startswith('#CHROM'): break

samps = pd.DataFrame(line.strip().split('\t')[9:], columns=['samps'])

#get sNMF Q file
dat = pd.read_csv(qfile, header=None, sep=' ') #FOLDER!!!!
dat = pd.concat([dat, samps], axis=1)

#get a pop file
pops = pd.read_csv('./pops.csv', header=0) 

dat = pd.merge(dat, pops, on='samps')
dat = dat.sort_values('population')
dat = dat.reset_index(drop=True)
pop = dat.population[0] #get first pop name from top of df

dat[['samps','population']].to_csv(f"{fldr}/{vcf.replace('.vcf','.structPlot')}_{qfile.split('.')[-2]}_Samples.csv", index=False)

im = np.full((length, width, 3), 245, dtype=np.int16) #full whitish rectangle
colDik = {'black': (0, 0, 0), 'red': (255, 0, 0), 'green': (0, 200, 0), 'cyan': (0,255,255), 'yellow': (255, 255, 0), 'purple': (128,0,128), 'orange': (255, 140, 0), 'silver': (192,192,192), 'teal': (0,128,128), 'blue': (200, 0, 0), 'dkred': (0, 0, 150), 'olive': (128,128,0), 'maroon': (128,0,0), 'grey': (128,128,128), 'white': (255,255,255), 'magenta': (255,0,255), 'lime': (0,255,0), 'navy': (0,0,128)}
keys = list(colDik.keys())[1:]

top = int(length * 0.5333) #top of panel
#Draw border - heights not flexible!! - could move to after where main plot is drawn
im[375: 400, w: width - w, :] = colDik.get('grey')
im[1995: 2025, w: width - w, :] = colDik.get('grey')
im[375: 2025, w - 25: w, :] = colDik.get('grey')
im[375: 2025, width - w: width - w + 35, :] = colDik.get('grey') #35 here not 25 as final stripes tend to overwrite border (little bit neater)

incrs = np.linspace(w, width - w, dat.shape[0]).astype(int) #boundaries btwn indvs
max_diff = np.diff(incrs).max() #get max diff btwn linspace gaps to make sure all space is covered (ie no gaps btwn bars) 
pop_limits = [w] #boundaries btwn pops + start posn (for calc midpoint placement of pop labels)

#draw main plot of Qvals
for i in dat.index:
    w = incrs[i]
    vals = np.array(dat.iloc[i,:])
    h = int(length * .13333) #h = top of panel
    cnt = 0
    for val in vals[:-2]: #DON'T cycle thru last two column values (these are samp and pop!!)
        im[h: h + int(val * top), w:w + max_diff, :] = colDik.get(keys[cnt]) #colour bar up to max len of dict
        h += int(val * top)
        cnt += 1
    if dat.population[i] != pop: #draw balck bar to seperate pops
        pop = dat.population[i]
        im[400: 2000, w - 10: w + 10, :] = colDik.get('black') #colour bar [not flexible!!!]
        pop_limits.append(w)


pop_limits.append(width - pop_limits[0]) #add last plot posn (calc midpoints of pop labels)
pops = dat.population.unique().tolist()

PIL_image = Image.fromarray(np.uint8(im)).convert('RGB')

#annotate pop labels
def rotateText(pop_limits, pops, angle, PIL_image): #d (depth) of text start point is not flexible
    d = int(.95* length) - (len(pops) * [150, 120][len(pops) > 6]) # contains a one line if then - a fudgey adjustment of pop label posn placement
    d = 2000 #!!!!!!!!!!!!!!!!!!!!!!
    font = ImageFont.truetype(sys.argv[3], int(640/len(pops))) # use a truetype font
    
    for pop in pops:
        fontColor = colDik.get('black') #colours reversed for pil image
        txt = Image.new('L', (500, 200), color=0) #mk text box
        imd = ImageDraw.Draw(txt)
        imd.text((0, 0), pop,  font=font, fill=255) #0=black
        w = txt.rotate(angle,  expand=1)
        x = int(((pop_limits[pops.index(pop) + 1] - pop_limits[pops.index(pop)]) / 2) + pop_limits[pops.index(pop)])
        PIL_image.paste(ImageOps.colorize(w, (0, 128, 0), fontColor), (x - int(width / 72), d),  w) #w, ?, colour, posn

    return PIL_image

PIL_image = rotateText(pop_limits, pops, 90, PIL_image)

PIL_image.save(f"{fldr}/{vcf.replace('.vcf','.structPlot')}_{qfile.split('.')[-2]}_nPops{len(pops)}.pdf", dpi=(600, 600)) #CHECK FOR sys!!!

