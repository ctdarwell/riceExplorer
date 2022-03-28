import glob, sys, numpy as np
import matplotlib.pyplot as plt

indir, outdir, job = sys.argv[1], sys.argv[2], sys.argv[3] #in out jobname

files = glob.glob1(indir, '*log')

#INDICES 1 & 2 BELOW!!!!
dat = np.array([[x.split('.')[1].replace('run', ''), x.split('.')[2].replace('K', '')] for x in files])
nruns = np.max(dat[:,0].astype(int))
kmax = np.max(dat[:,1].astype(int))
res = []

for k in range(1, kmax + 1):
    for n in range(1, nruns + 1):
        for line in open(f'{indir}/{job}.run{n}.K{k}.log'): #var=fldr
            if 'Cross-Entropy (masked data):' in line:
                res.append([k, line.split(' ')[-1].strip()])

res = np.array(res)
means = [np.mean(res[np.where(res[:, 0] == str(v)), 1].astype(float)) for v in range(1, kmax + 1)]
sem = [(np.std(res[np.where(res[:, 0] == str(v)), 1].astype(float)) / np.sqrt(len(means))) for v in range(1, kmax + 1)]

plt.errorbar(range(1, kmax + 1), means, yerr=sem, color='k', fmt="o")
plt.xlabel('No. of demes (k)')
plt.ylabel('Cross entropy (µ / SEM)')
plt.savefig(f"{outdir}/{job}.pdf")

#find K and choose file
cnt = 1
files = np.array(files)
for mu in means:
    if mu == min(means): k = f"K{cnt}.log"
    cnt += 1

opts = files[np.char.endswith(files, k)]
g = open('sNMF.txt','w')
g.write(np.random.choice(opts))
g.close()


