import pandas as pd, sys, numpy as np

prfx = sys.argv[1]
pops = pd.read_csv('pops.csv', header=0)
dik = dict(zip(pops.samps, pops.population))

g = open(f'{prfx}.allChroms.names.tre', 'w')
for line in open(f'{prfx}.allChroms.nex.tre'):
    for d in dik:
        line = line.replace(d, f"{d}_{dik.get(d)}")
    g.write(line)
g.close()

