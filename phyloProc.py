import pandas as pd, sys, numpy as np
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import concurrent.futures, itertools

'''NB RE METHODS
In order to statistically assess the number of likely clusters (K), Ward's method calculates a
merging cost (Δ; how much the sum of squares increases when two clusters are merged) at each step
when incrementally merging clusters from n (i.e. all values) to one (i.e. a single cluster).
Excessive merging costs are disqualified and the value of K+1 is taken to represent K
'''

fas = sys.argv[1] #fasta file
indir, outdir = sys.argv[2], sys.argv[3] 
max_workers = None #nProcesses - leave as default on HPC [ie None]; set as no. processors on your machine -1 for laptops/Desktops 
chunksize = 1 #break up of parallelisation chunks; large values may be faster for larger datasets

#start
taxa, seqs = [], []
for line in open(f"{indir}/{fas}"): 
    if line.startswith('>'): taxa.append('_'.join(line[1:-1].replace('.1','').split(' ')[:2]))
    #if line.startswith('>'): taxa.append(line[1:-1])
    else: seqs.append(line)

def main():
    mat = np.zeros([len(seqs), len(seqs)])
    dik = getDik() #k2p ts/tv probabilities for ambig codes
    
    #calc K2P dists in plll
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        pairs = [pair for pair in itertools.combinations(taxa, 2)]
        res = executor.map(k2p, pairs, itertools.repeat(dik), chunksize=chunksize)
    
    #extract results and format K2P matrix
    dists = [r for r in res]
    cntr = 0
    for pair in pairs:
        mat[taxa.index(pair[0]), taxa.index(pair[1])] = dists[cntr]
        mat[taxa.index(pair[1]), taxa.index(pair[0])] = dists[cntr]
        cntr += 1
    
    #further processing
    #write K2P distance matrix
    out = pd.DataFrame(mat)
    out.columns = taxa
    out.index = taxa
    out.to_csv(f"{outdir}/{fas.replace(fas.split('.')[-1], 'k2p.csv')}")
    
    #make NJ tree and plot
    Z = sch.linkage(mat, method='ward', optimal_ordering=True)
    plt.figure(figsize=(12,12))
    if len(taxa) > 80: dendogram = sch.dendrogram(Z, leaf_font_size=9, no_labels=True, orientation='left', distance_sort=True) #https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html
    else: dendogram = sch.dendrogram(Z, leaf_font_size=9, orientation='left', labels=np.array(taxa), distance_sort=True) #https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html
    plt.savefig(f"{outdir}/{fas.replace(fas.split('.')[-1], 'pdf')}")
    
    '''
    #write clusterings - see matplotlib colors C0 to C9:
    colDik = {'C0':'blue', 'C1':'orange', 'C2':'green', 'C3':'red', 'C4':'purple', 'C5':'brown', 'C6':'pink', 'C7':'gray', 'C8':'olive', 'C9':'cyan'}
    grps = pd.DataFrame([dendogram['leaves_color_list'], dendogram['ivl']]).T
    grps.columns = ['cluster','samps']
    cnt = 1
    if not np.isin([False], np.isin(grps.cluster.unique(), list(colDik.keys()))):
        for grp in grps.cluster.unique():
            grps.cluster[grps.cluster == grp] = colDik.get(grp)
            cnt += 1
    
    grps.samps = [taxa[s] for s in dendogram['leaves']]
    grps.to_csv(f"{outdir}/{fas.replace(fas.split('.')[-1], 'nj.clusters.csv')}", index=False)
    '''
    #make Newick
    tree = sch.to_tree(Z, False)
    nwk = get_newick(tree, tree.dist, taxa)
    
    #write Nexus tree
    writeTree(taxa, nwk)

def get_newick(node, parent_dist, leaf_names, newick='') -> str:
    """
    Convert sciply.cluster.hierarchy.to_tree()-output to Newick format.
    https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format

    :param node: output of sciply.cluster.hierarchy.to_tree()
    :param parent_dist: output of sciply.cluster.hierarchy.to_tree().dist
    :param leaf_names: list of leaf names
    :param newick: leave empty, this variable is used in recursion.
    :returns: tree in Newick format
    """
    
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parent_dist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parent_dist - node.dist, newick)
        else:
            newick = ");"
        newick = get_newick(node.get_left(), node.dist, leaf_names, newick=newick)
        newick = get_newick(node.get_right(), node.dist, leaf_names, newick=",%s" % (newick))
        newick = "(%s" % (newick)

        return newick

#make dict of ts-tv K2P distances inc probabilistic vals for ambig nuc codes
def getDik():
    dik = {}
    trans = [['A','G'], ['C','T']]
    ambs = {'A':'A','C':'C','T':'T','G':'G', 'R': 'AG', 'Y': 'CT', 'M': 'AC', 'W': 'AT', 'S': 'CG', 'K': 'GT'}

    for comb in itertools.combinations(ambs, 2):
        inst = [list(ambs.get(comb[0])), list(ambs.get(comb[1]))]
        l = []
        perms = list(itertools.product(inst[0], inst[1]))
        for ins in perms:
            if ins[0] == ins[1]: continue
            if sorted(ins) == trans[0] or sorted(ins) == trans[1]: l.append('ts')
            else: l.append('tv')
        dik[''.join(sorted(comb))] = [l.count('ts')/len(perms), l.count('tv')/len(perms)]
        
    return dik

#k2p ts/tv probabilities for ambig codes
def k2p(pair, dik):
    s1, s2 = seqs[taxa.index(pair[0])].upper(), seqs[taxa.index(pair[1])].upper()
    ts, tv = 0.0, 0.0
    for i in range(len(s1)):
        qw = ''.join(sorted([s1[i], s2[i]]))
        if qw in dik:
            ts += dik.get(qw)[0]
            tv += dik.get(qw)[1]

    # p = transition freq; q = transversion freq 
    p = ts / len(s1)
    q = tv / len(s1)
    
    #calc K2P - may give some neg vals for np.log    
    logv = (1 - (2 * p) - q) * np.sqrt(1 - (2 * q))
    if logv <= 0: kimura = 0.999999 #this is an ugly fudge [cf. R/ape generates an error] but should only apply for ridiculously distant sequences
    else: kimura = -0.5 * np.log(logv)
    
    return kimura


def writeTree(taxa, nwk):
    f = open(f"{outdir}/{fas.replace(fas.split('.')[-1], 'nex.tre')}", "w")
    f.write(f"#NEXUS\n\nBEGIN TAXA;\n\tDIMENSIONS NTAX = {len(taxa)};\n\tTAXLABELS\n")
    for tax in taxa:
        f.write(f"\t\t{tax}\n")
    f.write("\t;\nEND;\nBEGIN TREES;\n\tTRANSLATE\n")
    c = 1
    for tax in taxa:
        f.write(f"\t\t{c}\t{tax},\n")
        c += 1
    
    f.write("\t;\n\tTREE * UNTITLED = [&U]\n")
    
    for tax in taxa:
        nwk = nwk.replace(tax, str(taxa.index(tax) + 1))
    
    f.write(f"{nwk}\nEND;\n")
    f.close()


if __name__ == '__main__': main()

