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
10. Makes summary table of all accessions with extreme haplos and a counts X accession table (*pipeline_summary.py*)








