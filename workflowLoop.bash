#!/bin/bash

module load SAMtools/1.9-intel-2019b BCFtools/1.10.2-intel-2019b #MAFFT/7.453-iccifort-2019.5.281-with-extensions

fldr=chr03_GL 				#work folder
samps=irri102_GL_paths.csv		#samps in irri_nipp with trait data
data=irri102_GL_data.csv 		#corresponding trait data
column=GL 				#name of column in data file with trait data
genes=msu_chr03.csv		#List of genes being examined
accs=Rice_accession_full_list.csv	#List of samples being examined for presence/absence of id'd genes

#ADD ''' HERE IF NEED TO PART RUN WORFLOW

## *** START WORKFLOW ***
mkdir $fldr
rm $fldr/* #should be empty fldr as we write using '>>'

while read q; #gene list i.e. MSU genes
do echo GENE: $(echo $q | cut -d',' -f 2);

#DO bcftools
echo $'#'CHROM$'\t'POS$'\t'ID$'\t'REF$'\t'ALT$'\t'QUAL$'\t'FILTER$'\t'INFO$'\t'FORMAT$'\t'misc >> $fldr/$(echo $q | cut -d',' -f 2).vcf #write header line
while read p; #samples from samps
do echo $(echo $p | cut -d'/' -f 3 | cut -d'.' -f 1) >> $fldr/$(echo $q | cut -d',' -f 2).vcf #write sample name
#select all variants from a gene region across IRRI accs; --min-ac=1: "AC allele count in genotypes, for each ALT allele"
bcftools view --min-ac=1  $p $(echo $q | cut -d',' -f 1):$(echo $q | cut -d',' -f 4)-$(echo $q | cut -d',' -f 5) | sed -n '/AC=/,$p' >> $fldr/$(echo $q | cut -d',' -f 2).vcf #NB 'sed' looks for 1st 'AC=' and writes only this and following rows to file
done < $samps
done < $genes # msu_testDat.csv #msu_ver7_simple.csv #END read q e.g. msu_testDat.csv
echo 1. DONE: bcftools-IRRI

#rm repeat entries (i.e. from diff accs), format for snpEff & write allele X sample file: *.posn_samps.csv; sys.argv[2] denotes identifying prefix of sample name (to distinguish lines from rest of data)
python vcf4snpEffLoop2.py LOC*vcf $fldr IRIS #file search term, folder, sample prefix
echo 2. DONE: vcf4snpeff.py[1]

#run snpEff - id potentially influential SNPs
vcfs=$(ls $fldr/LOC*vcf)
for v in $vcfs
do java -Xmx4g -jar /tarafs/data/home/cdarwell/snpEff/snpEff.jar Oryza_sativa $v > $(echo $v | cut -d'.' -f 1).ann.vcf
#select only loci with significant predicted impact
grep 'MODERATE' $(echo $v | cut -d'.' -f 1).ann.vcf >> $(echo $v | cut -d'.' -f 1).ann.best
grep 'HIGH' $(echo $v | cut -d'.' -f 1).ann.vcf >> $(echo $v | cut -d'.' -f 1).ann.best
done # msu_testDat.csv #msu_ver7_simple.csv #END read q e.g. msu_testDat.csv
rm $fldr/LOC*.ann.vcf
echo 3. DONE: snpEff

#evaluate which IRRI haplotypes (based on snpEff influential SNPs), have extremes phenotypes 
python hapXphenoPredictorSNPsLoop2.py $column $data $fldr
echo 4. DONE: hapXphenoPredictorSNPs.py - checked haps!

#search RGD for ALT alleles in focal gene region: OUT: rgd_*.vcf
rgd=$(ls ../rgd_gvcf/*gz)
#while read q; #gene list i.e. MSU genes
loci=$(ls $fldr/*haps.csv)
for loc in $loci
do q=$(grep $(echo $loc | cut -d'/' -f 2 | cut -d'.' -f 1) $genes)
echo $'#'CHROM$'\t'POS$'\t'ID$'\t'REF$'\t'ALT$'\t'QUAL$'\t'FILTER$'\t'INFO$'\t'FORMAT$'\t'misc >> $fldr/rgd_$(echo $q | cut -d',' -f 2).vcf
for r in $rgd
do echo $(echo $r | cut -d'/' -f 3 | cut -d'.' -f 1) >> $fldr/rgd_$(echo $q | cut -d',' -f 2).vcf
bcftools view --min-ac=1 $r $(echo $q | cut -d',' -f 1):$(echo $q | cut -d',' -f 4)-$(echo $q | cut -d',' -f 5) | sed -n '/AC=/,$p' >> $fldr/rgd_$(echo $q | cut -d',' -f 2).vcf
done
done < $genes # msu_testDat.csv #msu_ver7_simple.csv #END read q e.g. msu_testDat.csv
echo 5. DONE: bcftools-RGD

#format bcftools file
python vcf4snpEffLoop2.py rgd_LOC*vcf $fldr W00 #file search term, folder, sample prefix
echo 6. DONE: vcf4snpeff.py[2]

#look for SNPs coding extreme phenotypes in IRRI data (e.g. hi or lo) that are present in landraces/RD lines that are not found in RGD varieties
python potential_rgd_sourcesLoop4.py $fldr $accs #OUT: *potential*donors.csv
echo 7. DONE: potential_rgd_sources.py

