#!/bin/bash
#SBATCH --cpus-per-task=24

#module load SAMtools/1.12-GCC-10.2.0 BCFtools/1.12-GCC-10.2.0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

fldr=GL	 		                  	## 1. work folder [to be created]
genes=msu_115testGenes.csv          ## 2. List of genes being examined - FORMAT FOLLWS MSU DOWNLOAD! (NB msu_ver7_simple.csv is a simplified version of downloaded annotated genes from MSU)
rgd=$(cat rgdSampPaths.csv)     	## 3. file containing paths to gvcf files of local panel being investigated (eg RGD panel)
samps=irri100_GL_paths.csv	     	## 4. paths to IRRI gvcf files - NO ESSENTIAL FORMAT
data=irri100_GL_data.csv 	      	## 5. corresponding IRRI trait data
column=GL 			              	## 6. name of column in IRRI trait data file
accs=Rice_accession_full_list.csv	## 7. List of examined samples (MUST have columns "acc" AND "type" - includes accession type: e.g. 'RGD improved line' vs. 'Landrace')
idd_genes=msu_GL.csv	           	## 8. simple list of known genes of interest for relevant phenotype (with NO header); put NA if know information is available 
all_msu=msu_ver7_simple.csv	    	## 9. reqd for figure output pgm for abridged test datasets so idd_genes are included
gene_id=LOC_Os                      ## 10. suffix of annotated genes (this is LOC from MSU); NB please READ GitHub notes
db_id=IRIS                          ## 11. suffices of database gvcf files (ALL files must start with same suffix) 
focal_id=W_OS_                      ## 12. suffices of investigated samples gvcf files (ALL files must start with same suffix) 
fig1="RD_variety,Landrace"          ## 13. The names of Non-focal variety (NFVs) types - IMPORTANT: use commas between types and replace spaces within type names with an underscore '_'
fig2="RGD_improved_line"            ## 14. The names of Elite cultivar (ECs) types - IMPORTANT: use commas between types and replace spaces within type names with an underscore '_'
db=Oryza_sativa                     ## 15. Name of snpEff database required
chrom=1                             ## 16. Chromosomes column number in [2]
loc_col=2                           ## 17. Annotated gene names column number in [2]
first=4                             ## 18. Column number of first base positions in [2]
end=5                               ## 18. Column number of end base positions in [2]

#ADD ''' HERE IF NEEDED TO PART RUN WORFLOW:

### *** START OF WORKFLOW ***
mkdir $fldr
rm $fldr/* #should be empty fldr as we write using '>>'

bcftask(){
while read p; #samples from samps
do echo $(echo $p | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1) >> $3/$(echo $1 | cut -d',' -f 2).vcf
bcftools view --min-ac=1  $p $(echo $1 | cut -d',' -f $4):$(echo $1 | cut -d',' -f $5)-$(echo $1 | cut -d',' -f $6) | sed -n '/AC=/,$p' >> $3/$(echo $1 | cut -d',' -f 2).vcf
done < $2
}

export -f bcftask
xargs -a $genes -rd'\n' -P20 -I'{}' bash -c 'bcftask "$1" "$2" "$3" "$4" "$5" "$6"' bash {} "$samps" "$fldr" "$chrom" "$first" "$end" #syntax note: {} represents $genes

echo 1. DONE: bcftools-IRRI
wait


#rm repeat entries (i.e. from diff accs), format for snpEff & write allele X sample file: *.posn_samps.csv; sys.argv[2] denotes identifying prefix of sample name (to distinguish lines from rest of data)
python vcf4snpEffLoop.py $gene_id*vcf $fldr $db_id #file search term, folder, sample prefix
echo 2. DONE: vcf4snpEffLoop.py-1
wait


#run snpEff - id potentially influential SNPs
snpefftask(){
java -Xmx4g -jar /fs/home/cliveterence.dar/snpEff/snpEff.jar $2 $1 > $(echo $1 | cut -d'.' -f 1).ann.vcf
#select only loci with significant predicted impact
grep 'MODERATE' $(echo $1 | cut -d'.' -f 1).ann.vcf >> $(echo $1 | cut -d'.' -f 1).ann.best
grep 'HIGH' $(echo $1 | cut -d'.' -f 1).ann.vcf >> $(echo $1 | cut -d'.' -f 1).ann.best
}

export -f snpefftask
find $fldr/ -name $gene_id*vcf | xargs -P20 -I'{}' -d'\n' bash -c 'snpefftask "$1" "$2"' bash {} $db #searches for all LOC*vcf files

rm $fldr/$gene_id*.ann.vcf
find $fldr/*best -size 0 -delete
echo 3. DONE: snpEff
wait

#evaluate which IRRI haplotypes (based on snpEff influential SNPs), have extremes phenotypes 
python hapXphenoPredictorSNPsLoop.py $column $data $fldr
echo 4. DONE: hapXphenoPredictorSNPs.py - checked haps!
wait


#search RGD for ALT alleles in focal gene region: OUT: rgd_*.vcf
bcftask2(){
q=$(grep $(echo $1 | cut -d'/' -f 2 | cut -d'.' -f 1) $3) #$3 = $genes
for r in $2
do echo $(echo $r | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1) >> $4/rgd_$(echo $q | cut -d',' -f $8).vcf
bcftools view --min-ac=1 $r $(echo $q | cut -d',' -f $5):$(echo $q | cut -d',' -f $6)-$(echo $q | cut -d',' -f $7) | sed -n '/AC=/,$p' >> $4/rgd_$(echo $q | cut -d',' -f $8).vcf
done
}

export -f bcftask2
find $fldr/ -name "*haps.csv" | xargs -P20 -I'{}' -d'\n' bash -c 'bcftask2 "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8"' bash {} "$rgd" "$genes" "$fldr" "$chrom" "$first" "$end" "$loc_col"

echo 5. DONE: bcftools-RGD
wait


#format bcftools file
python vcf4snpEffLoop.py rgd_$gene_id*vcf $fldr $focal_id #file search term, folder, sample prefix
echo 6. DONE: vcf4snpEffLoop.py-2


#look for SNPs coding extreme phenotypes in IRRI data (e.g. hi or lo) that are present in landraces/RD lines that are not found in RGD varieties
python potential_rgd_sourcesLoop.py $fldr $accs $fig1 #OUT: *potential*donors.csv & *assocSNPs.csv
echo 7. DONE: potential_rgd_sources.py


#make figure - NB pipelineGraphOut_bashLoop.py
python pipelineGraphOut.py $fldr $all_msu $idd_genes $gene_id $fig1 $fig2 $chrom $loc_col $first $end
echo 8. DONE: mk figure!

#summarise output files
python pipeline_summary.py $fldr $fig1
echo finito!

