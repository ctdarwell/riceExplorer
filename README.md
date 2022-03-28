# riceExplorer

If you need assistance: ctdarwell@gmail.com (Clive).
This is a Linux/Bash workflow to to identify haplotype reservoirs among a genomic rice panel (or other organism). The workflow depends on BCFtools, snpEff and Python3 installed to your Linux environement. You are required to download *.gvcf* files from a global database (e.g. rice; IRRI: http://iric.irri.org/resources/3000-genomes-project; aka 3KRG) representing global genetic and phenotypic diversity of the study organism. You require associated trait data for the files downloaded from the global database for your trait of interest. You then require a local genomic panel of accessions (*.gvcf*) that you want to mine for potentially beneficial SNPs/haplotypes that corellate with extreme phenotypes from the database. The pipeline initially identifies and clusters functional haplotypes of genome-wide gene regions from the global database and then calculates  phenotypic associations and evaluates whether they are above or below user-inputted quantile ranges. Subsequently, the pipeline searches your local genomic database to recognise if identified markers are found only among selected non-focal varieties (NFVs) versus elite cultivars (ECs). The basic premise is that most agricultural breeding programmes focus on ECs (which over time may may accumulate deleterious alleles) and may come to neglect potential contributions of other accessions held by the research body. It represents a fully automated bioinformatics pipeline that can easily integrate local breeding program sequence data with international database resources, without relying on any (costly, time-consuming) phenotypic experimental procedure. The intended result is that identified loci may serve as genetic donors to the development panel. **riceExplorer** can automatically conduct a full genome analysis and produces annotated graphical output of chromosomal maps, potential global diversity sources, and summary tables.

**Dependencies**

This resource comprises several Python3 programs that integrate via a ***Bash*** script with freely-available bioinformatics software in the public domain. In order to run it must be able to call both the **Bcftools** (https://samtools.github.io/bcftools/) and **snpEff** (http://pcingola.github.io/SnpEff/) program suites in a Linux environment. Additionally, you must set several parameters (data file names and other details) at the start of the Bash script so that the program may proceed.
PYTHON libraries required: pandas, numpy, PIL, scipy, seaborn, matplotlib, cv2 (normally: `pip install <library> --user`) on your HPC terminal
OPTIONALLY (though recommended) you must utilise R (normally found on HPC distributions). R libraries **LDheatmap** and **snpStats** must be installed.
OPTIONALLY you may download and install the software sNMF (http://membres-timc.imag.fr/Olivier.Francois/snmf/contact.htm) in you work folder for the full program to work although the pipeline will still perform phylogenetic, LD and site-frequency analyses without it.

**Preliminary setup**

The program must first evaluate the genomic variation encoded within downloaded files from a global (or other repository) database. For rice (*Oryza sativa*), the IRRI database hosts 3000+ whole-genome sequenced accessions that may be freely downloaded. A first step is to decide which samples to investigate that represent global diversity. For this we provide the auxilliary Python script *sampleSelector.py* (see below) that takes a table of accessions, their varietal type (for rice this includes *indica*, *japonica* etc), the geographic origin, and their associated phenotype for the trait of interest. *sampleSelector.py* then outputs a list of accessions that includes the largest and smallest phenotypic value for each varietal type for each global region. From this, representative accessions of the global panel should include the widest spread of phenotypic and genotypic diversity. After this process the full program is ready to run. It should be noted that no phenotypic data is necessary for the local panel (i.e. accessions under investigation by a breeding program team). Additionally, while SNPs from both panels should be called against the same reference genome (e.g. for rice, Nipponbare), there should also be consistent naming of chromosomes across panels. For example, a reference genome may have chromosomes referred to as either *>Chr01* or *>1* (i.e. in a fasta formatted file) - care should be taken that subsequent *GVCF* files used by **riceExplorer** have consistent chromosome calling.

**MAIN WORK-FLOW**
1. PYTHON evaluates global database to select accessions of interest (*sampleSelector.py*)
2. BCFTOOLS records all SNPs from selected global database accessions 
3. PYTHON formats files for SNPEFF and records which SNP are associated with which global database accession (*vcf4snpEffLoop.py*)
4. SNPEFF evaluates functional importance of SNPs (inc BASH filtering)
5. PYTHON outputs significant SNPs (of extreme values) with associated global database accessions (*hapXphenoPredictorSNPsLoop.py*)
6. BCFTOOLS records all SNPs from selected local panel accessions (at loci that have identified extreme value SNPs)
7. PYTHON formats files (as for SNPEFF) (*vcf4snpEffLoop.py*)
8. PYTHON outputs SNPs coding extreme phenotypes in global data (e.g. hi or lo) that are present ONLY in NFVs and also outputs all identified extreme haplotypes (*potential_rgd_sourcesLoop.py*)
9. PYTHON makes graphs of selected NFVs vs. ECs (*pipelineGraphOut.py*)
10. PYTHON makes summary table of all accessions with extreme haplos and a counts X accession table (*pipeline_summary.py*)

ADDITIONAL: it should be noted that running steps 2-5 on a local panel (i.e. of research samples) can highlight SNPs haplotypes associating with extreme phenotypes without comparing them against a 2nd panel (this may also be useful)

**Variable values required for *riceExplorer_vars***. While care must be taken for the program to run, this is all straightforward, non-technical information.

1. VARIABLE: `fldr` - input a folder name that will be created by the program to output results and conduct analyses
2. VARIABLE: `genes` - name of file containing names of all genes (gene regions) to be investigated including columns indicating chromosome numbers, gene annotation names, start and end base positions to be examined. **NB** remove column headers [see "msu_115testGenes.csv" in *dataFiles.zip*]
3. VARIABLE: `rgd` - name of file containing paths of all local panel *.gvcf* files [see "rgdSampPaths.csv" in *dataFiles.zip*]
4. VARIABLE: `samps` - name of file containing paths of all global panel *.gvcf* files [see "irri200_GL_paths.csv" in *dataFiles.zip*]
5. VARIABLE: `data` - name of file containing trait data for each *.gvcf* file sample (MUST have column name "acc" & "trait" for accession and phenotypic value) [see "irri100_GL_data.csv" in *dataFiles.zip*]
6. VARIABLE: `accs` - name of file of list containing examined samples (MUST have columns "acc" AND "type" - includes accession type: e.g. 'Landrace', 'Dept variety') [see "Rice_accession_full_list.csv" in *dataFiles.zip*]
7. VARIABLE: `idd_genes` - name of file of list of genes known to have functional impact on relevant phenotype; NO header, NB put NA if know information is available [see "msu_GL.csv" in *dataFiles.zip*]
8. VARIABLE: `all_msu` - name of file of same format as `genes` - may often be same file; required in conjunction with `idd_genes` to annotate those genes on graphical output [see "msu_115testGenes.csv" in *dataFiles.zip*]
9. VARIABLE: `gene_id` - this is a search term used by the program to identify files assocaited with each gene; e.g. the Michigan State University annotation of the Nipponbare rice genome uses notation as: *LOC_Os01g02700*. It is intended that gene notation follows a similar scheme. Here there is a repeatable suffix ("*LOC_Os*") preceding every chromososme number which are then followed by a letter *g* that delimits the number of the gene; for such MSU genes we use: `gene_id=LOC_Os`
10. VARIABLE: `db_id` - a search term repressenting the suffix given to all global database *.gvcf* files (e.g. `db_id=IRIS` for IRRI samples)
11. VARIABLE: `focal` - a search term repressenting the suffix given to all local accession *.gvcf* files (e.g., all samples begin with "W00" for our data)
12. VARIABLE: `fig1` - the names of non-focal variety (NFVs) types (as listed under variable in [6]). This **MUST** use quotes and feature commas between words **WITHIN** names and underscores **BETWEEN** names (e.g. `fig1="RD_variety,Landrace"`) - i.e. there **MUST** be no spaces.
13. VARIABLE: `fig2` - the names of elite cultivar (ECs) variety types. Same as [12]: (e.g. `fig1="RGD_improved_line"`) - again, quotes, commas, underscores and **NO** spaces
14. VARIABLE: `db` - database called in `snpEff`. e.g. `db=Oryza_sativa`
15. VARIABLE: `chrom` - Chromosomes column number in [2] (e.g. `chrom=1`)
16. VARIABLE: `loc_col` - Annotated gene names column number in [2] (e.g. `loc_col=2`)
17. VARIABLE: `first` - Column number of first base positions in [2] (e.g. `first=4`)
18. VARIABLE: `end` - Column number of last base positions in [2] (e.g. `end=5`)
19. VARIABLE: `nChroms` - No. of chromosomes for the organism (e.g. `nChroms=12` 
20. VARIABLE: `thresh` - NB Two input options (e.g. `thresh=0.866`: 1. A percentile value - thresh=0.975 gives 5% confidence intervals; thresh=0.866 gives percentiles at 1.5 standard deviations; 2. lower and upper absolute values SEPERATED by a comma - e.g. a lower value of 7 and upper value 10 should be written as: thresh=7,10)  
21. VARIABLE: `font1` - Header font - on HPC look in directory '/usr/share/fonts/' for available fonts (e.g. `font1=FreeSansBold.ttf`)
22. VARIABLE: `font2` - minor font - on HPC look in directory '/usr/share/fonts/' for available fonts (e.g. `font2=FreeSans.ttf`)
23. VARIABLE: `LDtypes` - search term in pops files; searches for samples names used in calculating LD (e.g. `LDtypes=IND`) - searches for indica samples in our data (IND_IRRI or IND_RGD) - see below in `riceExplorer_aux.bash`


To run the workflow you should type: `./riceExplorer.bash` in a Linux terminal. Or, preferably you should run it through your system job loader (e.g., slurm) - the full rice analysis will take 2-3 days on a HPC. NB line 5 in the bash script (*riceExplorer.bash*) will need altering to account for the versions of BCFtools you are using. 

**Auxiliary analyses**
After running `riceExplorer.bash`, you may then run `riceExplorer_aux.bash`. It also uses `riceExplorer_VARS` so there is no need to reconfigure any data files. 
It will perform numerous diversity analyses (see paper)
NB you **must** include a file called **pops.csv** with header columns: **samps** and **population**. The **population** column should consider variable #23 `LDtypes`. For example, in our data we have indica samples from two panels. They are coded as "IND_IRRI" or "IND_RGD" - we can use `LDtypes=IND` to search for both (Please ensure the term does not feature in any sample names)

**AUXILIARY WORK-FLOW**
1. BASH evaluates high quality/good coverage SNPs
2. BCFTOOLS calls SNPs, creates VCF files
3. PYTHON tidies up VCF files
4. PYTHON converts VCFs to Fasta format
5. PYTHON performs phylogenetic analyses and creates tree file
6. PYTHON creates SNAPP and DIYABC input files
7. PYTHON generates site-frequency spectrum graphical output
8. R calculates LD relationships
9. PYTHON generates annotated LD figures
10. sNMF perform population genomic analyses


**Note on use of *sampleSelector.py***

For example, IRRI 3kRG has approximately 2000 accessions with linked grain length data. Use *sampleSelector.py* to reduce these samples to represent most phenotypic diversity from each accession type (e.g. for rice: *indica*, *japonica* etc) in each geographic region. The CSV file must contain the column names: *acc* [accession names]; *type* [accession types, e.g. *indica*, *japonica* etc]; region [e.g. Africa, China etc]; *trait* [phenotypic measurement]. A column named *country* may also be included in order to exclude any particular countries [e.g. the accession panel's native country].

Please use as: `python sampleSelector.py arg1` - where arg1 is the data file.

Optionally, use as `python sampleSelector.py arg1 arg2 arg3` - where arg2 is countries to be excluded (using commas, underscores and no spaces - e.g. Thailand,South_Africa), and arg3 is a separate file with accessions (column name: *acc*) to be excluded. If no countries are to be excluded replace arg2 with NA or XXX (or some other meaningless placeholder)

IRRI 3KRG has accession types *indica*, *japonica*, *aus*, *admix* and *aro*. However, *indica* is subdivided into *ind1A*, *ind1B* ect and *japonica* into *trop* and *temp* - you may wish to use such subdivisions under the main data file *type* column.



**Note on stand-alone use of *pipelineGraphOut.py***

This program can be used if you wish to make a seprate figure representing a specific sample - i.e. to indicate the potentially useful marker positions associated with a particular accession (e.g. Figure 6 in the paper). First, you must compile an extra CSV file. No need for a header. Each row should represent genes with useful haplotypes found on that sample, e.g. LOC_Os07g19530, LOC_Os07g25410 etc. Then, immediately underneath you may place any other gene names that are known to be significant on the chromosome being examined (same format of gene names). The program takes 14 arguments. These are [1] path to work folder; [2] List of genes being examined (equivalent to `genes`); [3] newly compiled filename; [4] gene suffix search term (e.g. "LOC_Os', equivalent to `gene_id`); [5] names  of non-focal variety (NFVs) types (equivalent to `fig1`); [6] names  of elite cultivar types (equivalent to `fig2`); [7] column number in list of genes for chromosomes; [8] column number in list of genes for gene names; [9] column number in list of genes for start position of gene; [10] column number in list of genes for end position of gene; [11] number of chromosomes for the organism; [12] name of sample (or other value to print on figure); [13] chromosome number being examined; [14] row number in newly compiled file where first known gene appears (e.g., if the first six genes in the file represent genes with high value haplotypes on that sample, and the next genes are known genes important for the phenotype being examined, you should add "7")

To run the program you should type: `python pipelineGraphOut.py arg1 arg2 arg3 arg4 arg5 arg6 arg7 arg8 arg9 arg10 arg11 arg12 arg13 arg14` in a Python terminal
