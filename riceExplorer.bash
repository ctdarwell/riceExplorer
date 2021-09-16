#!/bin/bash
#SBATCH --cpus-per-task=24

#The next line should be altered to specifically identify that versions of BCFtools you are using
module load BCFtools

source riceExplorer_vars #input parameters file
rgd=$(cat $paths)

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
python hapXphenoPredictorSNPsLoop.py $data $fldr $percentiles
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
python pipelineGraphOut.py $fldr $all_msu $idd_genes $gene_id $fig1 $fig2 $chrom $loc_col $first $end $nChroms
echo 8. DONE: mk figure!

#summarise output files
python pipeline_summary.py $fldr $fig1
echo finito!

