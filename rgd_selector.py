import sys

fldr = sys.argv[1]
prfx = sys.argv[2]

f = open('tomove.txt')
files = f.read()
f.close()

g = open('rgds.txt')
files2 = g.read()
g.close()

tomove = files.strip().split(' ')
rgds = files2.strip().split(' ')

s = ""
for vcf in tomove:
    if vcf.replace(fldr, f"{fldr}{prfx}") not in rgds:
        s += f"{vcf.replace(fldr, '').replace('.vcf', '')} "

g = open('rgd_toAugment.txt', 'w')
vcfs = s[:-1].split(' ')
for vcf in vcfs:
    g.write(vcf + '\n')
g.close()

print(vcfs)

