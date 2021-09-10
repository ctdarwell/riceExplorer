# riceExplorer

This is a Linux/Bash workflow to to identify haplotype reservoirs among a local genomic rice panel (or other organism). The workflow depends on BCFtools, snpEff and Python3 installed to your Linux environement. You are required to download *.gvcf* files from a global database (e.g. rice; IRRI: http://iric.irri.org/resources/3000-genomes-project) representing global genetic and phenotypic diversity of the study organism. You require associated trait data for the files downloaded from the global database for your trait of interest. You then require a local genomic panel of accessions (*.gvcf*) that you want to mine for potentially beneficial SNPs/haplotypes that corellate with extreme phenotypes from the database. The pipeline initially identifies and clusters functional haplotypes of genome-wide gene regions from the global database and then calculates mean phenotypic associations and evaluates whether they are above or below statistical confidence interval criteria. Subsequently, the pipeline searches your local genomic database to recognise if identified markers are found only among selected non-focal varieties (NFVs) versus elite cultivars (ECs). The basic premise is that most agricultural breeding programmes focus on ECs (which over time may may accumulate deleterious alleles) and may come to neglect potential contributions of other accessions held by the research body. It represents a fully automated bioinformatics pipeline that can easily integrate local breeding program sequence data with international database resources, without relying on any (costly, time-consuming) phenotypic experimental procedure. The intended result is that identified loci may serve as genetic donors to the development panel. **riceExplorer** can automatically conduct a full genome analysis and produces annotated graphical output of chromosomal maps, potential global diversity sources, and summary tables.

**Dependencies**

This resource comprises several Python3 programs that integrate via a ***Bash*** script with freely-available bioinformatics software in the public domain. In order to run it must be able to call both the **Bcftools** (https://samtools.github.io/bcftools/) and **snpEff** (http://pcingola.github.io/SnpEff/) program suites in a Linux environment. Additionally, you must set several parameters (data file names and other details) at the start of the Bash script so that the program may proceed.

**Preliminary setup**

The program must first evaluate the genomic variation encoded within downloaded files from a global (or other repository) database. For rice (*Oryza sativa*), the IRRI database hosts 3000+ whole-genome sequenced accessions that may be freely downloaded. A first step is to decide which samples to investigate that offer represent global diversity. For this we provide the auxilliary Python script *sampleSelector.py* that takes a table of accessions, their varietal type (for rice this includes *indica*, *japonica* etc), the geographic origin, and their associated phenotype for the trait of interest. *sampleSelector.py* then outputs a list of accessions that includes the largest and smallest phenotypic value for each varietal type for each global region. From this, representative accessions of the global panel should include the widest spread of phenotypic and genotypic diversity. After this process the full program is ready to run. It should be noted that no phenotypic data is necessary for the local panel (i.e. accessions under investigation by a breeding program team).

 **MAIN WORK-FLOW**
1. PYTHON evaluates global database to select accessions of interest (*sampleSelector.py*)
2. BCFTOOLS records all SNPs from selected IRRI accessions 
3. PYTHON formats files for SNPEFF and records which SNP are associated with which IRRI accession (*vcf4snpEffLoop.py*)
4. SNPEFF evaluates functional importance of SNPs (inc BASH filtering)
5. PYTHON outputs significant SNPs (of extreme values) with associated IRRI accessions (*hapXphenoPredictorSNPsLoop.py*)
6. BCFTOOLS records all SNPs from selected RGD accessions
7. PYTHON formats files (as for SNPEFF) (*vcf4snpEffLoop.py*)
8. PYTHON outputs SNPs coding extreme phenotypes in global data (e.g. hi or lo) that are present ONLY in NFVs and also outputs all identified extreme haplotypes (*potential_rgd_sourcesLoop.py*)
9. PYTHON makes graphs of selected NFVs vs. ECs (*pipelineGraphOut.py*)
10. PYTHON makes summary table of all accessions with extreme haplos and a counts X accession table (*pipeline_summary.py*)


**Variable values required for *pllWorkflowLoop.bash***. While care must be taken for the program to run, this is all straightforward, non-technical information.

1. VARIABLE: `fldr` - input a folder name that will be created by the program to output results and conduct analyses
2. VARIABLE: `genes` - filename containing names of all genes (gene regions) to be investigated including columns indicating chromosome numbers, gene annotation names, start and end base positions to be examined. Following Michigan State University format: [WEBPAGES]
3. VARIABLE: `rgd` - filename containing paths of all local panel *.gvcf* files
4. VARIABLE: `samps` - filename containing paths of all global panel *.gvcf* files
5. VARIABLE: `data` - filename containing trait data for each panel *.gvcf* file sample (MUST have column name *acc* & *trait* for accession and phenotypic value)
6. VARIABLE: `accs` - filename of list containing examined samples (MUST have columns "acc" AND "type" - includes accession type: e.g. 'RGD improved line' vs. 'Landrace')
7. VARIABLE: `idd_genes` - filename of list of genes known to have functional impact on relevant phenotype; NO header, NB put NA if know information is available 
8. VARIABLE: `all_msu` - filename of same format as `genes` - may often be same file; required in conjunction with `idd_genes` to annotate those genes on graphical output
9. VARIABLE: `gene_id` - this is a search term used by the program to identify files assocaited with each gene; e.g. the Michigan State University annotation of the Nipponbare rice genome uses [WHAT?] notation for genes, such as: *LOC_Os01g02700*. It is intended that gene notation follows a similar scheme. Here there is a repeatable suffix ("*LOC_Os*") preceding every chromososme number which are then followed by a letter *g* that delimits the number of the gene; for such MSU genes we use: `gene_id=LOC_Os`
10. VARIABLE: `db_id` - a search term repressenting the suffix given to all global database *.gvcf* files (e.g. `db_id=IRIS`)
11. VARIABLE: `focal` - a search term repressenting the suffix given to all local accession *.gvcf* files
12. VARIABLE: `fig1` - the names of non-focal variety (NFVs) types (as listed under variable in [6]). This **MUST** be added with commas between words **WITHIN** names and underscores **BETWEEN** names (e.g. `fig1="RD_variety,Landrace"`) - i.e. there **MUST** be no spaces.
13. VARIABLE: `fig2` - the names of elite cultivar (ECs) variety types. Same as [12]: (e.g. `fig1="RGD_improved_line"`) - again, commas, underscores and **NO** spaces
14. VARIABLE: `db` - databse called in `snpEff`. e.g. `db=Oryza_sativa`
15. VARIABLE: `chrom` - Chromosomes column number in [2] (e.g. `chrom=1`)
16. VARIABLE: `loc_col` - Annotated gene names column number in [2] (e.g. `loc_col=2`)
17. VARIABLE: `first` - Column number of first base positions in [2] (e.g. `first=4`)
18. VARIABLE: `end` - Column number of last base positions in [2] (e.g. `first=5`)

