#!/bin/bash
#SBATCH --cpus-per-task=24

source riceExplorer_vars
mkdir $fldr/tmp_fldr

rgd=$(cat $paths)
n=$(wc -l $paths | cut -d' ' -f 1)
min_data=$(($n * 8 / 10))
maf=$(($n * 5 / 100))

#find hiQual SNPs
bcftask3(){
	#q=$(grep $(echo $1 | cut -d',' -f 2) $3) #$3 now = $chrom_no
	q=$1
	file=$4/tmp_fldr/Chr$3_$(echo $q | cut -d',' -f $8).vcf
	
	for r in $2
	do bcftools view --max-ac 2 $r $(echo $q | cut -d',' -f $5):$(echo $q | cut -d',' -f $6)-$(echo $q | cut -d',' -f $7) | grep -v '#' | grep -v 'LowQual' >> $file
	done
	
	#Find good SNPs
	snps=$(cat $file | cut -d$'\t' -f 2) #list all SNPs
	list=$(printf '%s\n' "${snps[@]}" | sort | uniq -c | sort -k1,1nr -k2) #sort SNPs

	i=2
	j=1
	
	while [ $(echo $list | cut -d' ' -f $j) -gt "$n" ]
	do
	snp=$(echo $list | cut -d' ' -f $i)
	cnt=$(grep $snp $file | wc -l | cut -d' ' -f 1) #count SNP occurrences
	bad=$(grep $snp $file | grep 'LowQual' | wc -l)
	good=$(( $cnt - $bad ))
	
	if [ $good < "$n" ]
	then
		continue
	fi
	
	if [ $cnt -lt "$n" ]
	then
		continue
	fi
	
	#count variants 
	cnt1=$(grep $snp $file | grep '0/0' | wc -l)
	cnt2=$(grep $snp $file | grep '1/1' | wc -l)

	if [ $cnt1 -gt "$maf" ] && [ $cnt2 -gt "$maf" ]
	then
		echo $(echo $3 | sed 's/^0*//'):$snp >> Chr$3.good_snps.txt
		return
	fi
	(( i += 2 ))
	(( j += 2 ))
	
	done
}

export -f bcftask3

for chrom_no in $( eval echo {01..$nChroms} )
do
touch Chr$chrom_no.good_snps.txt

grep $gene_id$chrom_no $genes >> Chr$chrom_no.msu.csv #mk file from main gene list for locus only
chrom_genes=Chr$chrom_no.msu.csv

xargs -a $chrom_genes -rd'\n' -P24 -I'{}' bash -c 'bcftask3 "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8"' bash {} "$rgd" "$chrom_no" "$fldr" "$chrom" "$first" "$end" "$loc_col"

done

echo btask3b_finished >> report.txt


#ORDER hiQual SNP FILES
for chrom_no in $( eval echo {01..$nChroms} )
do
file=Chr$chrom_no.good_snps.txt
snps=$(cat $file | cut -d':' -f 2)
pad=$(for snp in $snps; do printf "%08d" $snp; echo " "; done)
s=$(echo ${pad[@]} | awk 'BEGIN{RS=" ";} {print $1}' | sort)

for snp in $s; do echo $(echo $chrom_no | sed 's/^0*//'):$(echo $snp | sed 's/^0*//') >> Chr$chrom_no.neworder.txt; done

done

#generate VCFs
rm cmb_paths.txt
cat $paths >> cmb_paths.txt
cat $samps >> cmb_paths.txt
cmb_paths=cmb_paths.txt
mkdir $fldr/GOOD_SNPS_HQ

bcftask4(){
prfx=$(echo $1 | cut -d'.' -f 2)
while read p; #samples from samps
do bcftools view $p $(cat $1) | grep -v 'LowQual' | bgzip >  $fldr/GOOD_SNPS_HQ$prfx.$(echo $p | rev | cut -d'/' -f 1 | rev | cut -d'.' -f 1).vcf.gz
done < $2
}

export -f bcftask4
find ./ -maxdepth 1 -name "Chr*.neworder.txt" | xargs -d'\n' -P20 -I'{}' bash -c 'bcftask4 "$1" "$2"' bash {} "$cmb_paths" #syntax note: {} represents $genes

for chrom_no in $( eval echo {01..$nChroms} )
do
files=$(ls $fldr/GOOD_SNPS_HQ/Chr$chrom_no*gz)
for file in $files; do tabix $file; done

bcftools merge --merge all $files > Chr$chrom_no.Mge.vcf

done

echo bcftask4_finished >> report.txt


#PERFORM AUXILIARY ANALYSES
mkdir phyProc

for chrom_no in $( eval echo {01..$nChroms} ); do python vcfPrcsr.py $chrom_no; done

mv *Mge.vcf phyProc/
mv *neworder.txt phyProc/
rm tmp*vcf


python vcfCombnr.py $fldr

mv *faoCmbn.vcf phyProc/


python vcf2fasta.py $fldr

python phyloProc.py $fldr.allChroms.fas . . 

python treeRenamer.py $fldr


python vcf2diyabc.py $fldr.allChroms.vcf . . pops.csv #CAN BE rx_vars variable name

python vcf2sfs.py $fldr.allChroms.vcf . .

python vcf2snapp.py $fldr.allChroms.vcf . .

mkdir figures
mv *pdf figures/

python vcfStripper.py $LDtypes #CHECK WHAT IS ON SERVER!!!!!!!!!!!!!!!!!!!

mv *faoLD.vcf phyProc/


for chrom_no in $( eval echo {01..$nChroms} ); do Rscript LDdecay.R $chrom_no; done

mv *$LDtypes*LD.vcf phyProc/


for chrom_no in $( eval echo {01..$nChroms} )
do for gnz in _hignz1.csv _lognz1.csv _hignz2.csv _lognz2.csv
	do python LDheatmapMaker_cbar.v3.0beta.py Chr$chrom_no.$LDtypes.LD.csv $fldr $gene_id$chrom_no$gnz $genes $gene_id $idd_genes $loc_col $first $end $font2
	done
done


mv *pdf figures/
mv *png figures/


python LD_LGH_tabber.py $fldr . $gene_id $genes $idd_genes

mv *$LDtypes.LD.csv phyProc/


./sNMF_CL_v1.2/bin/vcf2geno $fldr.allChroms.vcf sNMF_CL_v1.2/$fldr.geno

for i in 1 2 3 4 5 6 7 8 9 10 #n reps
do for j in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 #max no. of demes
do ./sNMF_CL_v1.2/bin/sNMF -x ./sNMF_CL_v1.2/$fldr.geno -K $j -q ./sNMF_CL_v1.2/$fldr.run$i.K$j.Q -g ./sNMF_CL_v1.2/$fldr.run$i.K$j.G -s $RANDOM -c > ./sNMF_CL_v1.2/$fldr.run$i.K$j.log
done
done


python sNMF_Xenter.py ./sNMF_CL_v1.2 . $fldr #indir outdir job_name


snmf=$(cat sNMF.txt)
python sNMF_annotater.py $snmf $fldr $font2
mv sNMF.txt phyProc/
mv *pdf figures/
mv $fldr.* phyProc/

mkdir sNMF_CL_v1.2/$fldr
mv sNMF_CL_v1.2/*.run* sNMF_CL_v1.2/$fldr/

echo phylo_finished >> report.txt


