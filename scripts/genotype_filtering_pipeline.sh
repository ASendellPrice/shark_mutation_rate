#!/bin/bash
#SBATCH -A snic2021-5-8
#SBATCH -p core -N 1
#SBATCH -t 1-00:00:00
#SBATCH -J MergeVCFs

#Activate conda environment and load required modules
conda activate GenomicsGeneral
module load bioinfo-tools bcftools BEDTools/2.29.2 samtools/1.5

#Create list of input VCFs
ls ../joint_genotyping/block_*[0-9].vcf.gz > vcf.list

#Merge VCF files, sort and index output
bcftools concat -f vcf.list -o ../joint_genotyping/Merged.vcf.gz -O z
bcftools sort -o ../joint_genotyping/Merged.sorted.vcf.gz -O z ../joint_genotyping/Merged.vcf.gz
bcftools index ../joint_genotyping/Merged.sorted.vcf.gz
#Count = 

#Filter merged VCF keeping only biallelic SNPs from non-repeat, highly mappable regions
bedtools intersect -a ../joint_genotyping/Merged.sorted.vcf.gz \
-b ../HiglyMappable_Unmasked_ranges.txt -header \
| bcftools view --exclude-types indels \
| bcftools view --max-alleles 2 \
| bgzip > Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.vcf.gz
#Count = 

#Convert to "geno" format for use later
python /proj/snic2020-2-19/private/shark/users/ash/BIN/genomics_general/VCF_processing/parseVCF.py \
-i Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.vcf.gz > Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.genotypes.txt

#Set sample genotype to missing if sample genotype quality < 20, as below Q20
#accuracy decreases rapidly.
#Illustrated nicely here: https://gatk.broadinstitute.org/hc/en-us/articles/360035531872
bcftools query -l Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.vcf.gz\
> all.samples
zcat Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.vcf.gz \
| bcftools filter -S . -e 'FMT/GQ[@all.samples] < 20' \
| bgzip > Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.vcf.gz
#Count = 

#Remove sites where parental genotypes are missing as these aren't informative.
zcat Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.vcf.gz \
| bcftools view -e 'GT[9]="mis" || GT[10]="mis"' \
| bgzip > Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing.vcf.gz
#Count = 

#Set sample genotype to missing if allele balance is:
#< 0.2 for heterozygote calls
#< 0.9 for homozygote calls 
#This is implemented using the jvarkit tool "vcffilterjdk"
source ../scripts/filter_het_homo_by_AB.sh \
Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing.vcf.gz \
Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing_AB-Filtered.vcf.gz

#Create a subset of snps where parents are homozygous for different alleles and all
#offspring are heterozygous and all sample genotypes are non-missing. These "known" heterozygous sites will
#be used as the basis for defining custom SNP filtering criteria.
bcftools query -l Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing_AB-Filtered.vcf.gz \
| head -n 9 > offspring.samples
zcat Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing_AB-Filtered.vcf.gz \
| bcftools view -i 'GT[9]="ref" & GT[9]="hom"' \
| bcftools view -i 'GT[10]="alt" & GT[10]="hom"' \
| bcftools view -i 'GT[@offspring.samples]="het"' \
| bcftools view -e 'GT[]="mis"' \
| bgzip > temp1.vcf.gz
zcat Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing_AB-Filtered.vcf.gz \
| bcftools view -i 'GT[9]="alt" & GT[9]="hom"' \
| bcftools view -i 'GT[10]="ref" & GT[10]="hom"' \
| bcftools view -i 'GT[@offspring.samples]="het"' \
| bcftools view -e 'GT[]="mis"' \
| bgzip > temp2.vcf.gz
bcftools index temp1.vcf.gz
bcftools index temp2.vcf.gz
bcftools concat --allow-overlaps temp1.vcf.gz temp2.vcf.gz \
| bgzip > Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing.offspringHet.vcf.gz
rm temp1.vcf.gz
rm temp2.vcf.gz
#Count = 374540 

#Based on the distribution of SNP quality annotations in the above VCF, calculate
#the filtering criteria (lower = 5th percentile, upper = 95th percentile).
#Filtering will be applied to the following site annotations:
#total read depth: DP
#mapping quality: MQ (lower bound only)
#mapping quality rank sum: MQRankSum
#base quality rank sum: BaseQRankSum
#read position rank sum: ReadPosRankSum
#quality by depth: QD (lower bound only)
#Filtering will then be applied to the following individual annotations:
#genotype quality: GQ (lower bound only)
#sample read depth: DP 
Rscript ../generate.custom.filters.R \
Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing.offspringHet.vcf.gz \
Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing.vcf.gz \
Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing.customFiltered.vcf.gz \
11

#Apply custom filtering to informative variants
source apply.custom.filters.sh
#Count = 

#From the custom filtered dataset, extract positons where both parents are homozygous
#for the REF allele and at least 1 offsping is heterozygous
zcat Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing.customFiltered_noSiteDP.vcf.gz \
| bcftools view -e 'GT[9]="mis" || GT[10]="mis"' \
| bcftools view -i 'GT[9]="ref" & GT[9]="hom"' \
| bcftools view -i 'GT[10]="ref" & GT[10]="hom"' \
| bcftools view -i 'GT[0]="het" || GT[1]="het" || GT[2]="het" || GT[3]="het" || GT[4]="het" || GT[5]="het" || GT[6]="het" || GT[7]="het" || GT[8]="het"' \
| bgzip > Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing.customFiltered.parentsHomoRef_noSiteDP.vcf.gz

#From the custom filtered dataset, extract positons where both parents are homozygous
#for the ALT allele and at least 1 offsping is heterozygous
zcat Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing.customFiltered_noSiteDP.vcf.gz \
| bcftools view -e 'GT[9]="mis" || GT[10]="mis"' \
| bcftools view -i 'GT[9]="alt" & GT[9]="hom"' \
| bcftools view -i 'GT[10]="alt" & GT[10]="hom"' \
| bcftools view -i 'GT[0]="het" || GT[1]="het" || GT[2]="het" || GT[3]="het" || GT[4]="het" || GT[5]="het" || GT[6]="het" || GT[7]="het" || GT[8]="het"' \
| bgzip > Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing.customFiltered.parentsHomoAlt_noSiteDP.vcf.gz

#Combine the last two VCFs into a single file, these represent high confidence putative mutations
bcftools index Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing.customFiltered.parentsHomoRef_noSiteDP.vcf.gz
bcftools index Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing.customFiltered.parentsHomoAlt_noSiteDP.vcf.gz
bcftools concat --allow-overlaps \
Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing.customFiltered.parentsHomoRef_noSiteDP.vcf.gz \
Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.belowGQ20setmissing.noParentalMissing.customFiltered.parentsHomoAlt_noSiteDP.vcf.gz \
| bgzip > Putative_Mutations.vcf.gz
#Count = 

#For putative mutations, set sample genotype to missing if allele balance is:
#< 0.2 for heterozygote calls
#< 0.9 for homozygote calls 
#This is implemented using the jvarkit tool "vcffilterjdk"
source ../scripts/filter_het_homo_by_AB.sh \
Putative_Mutations.vcf.gz \
Putative_Mutations.AB-Filtered.vcf.gz

#Convert putative mutations to "geno" format
python /proj/snic2020-2-19/private/shark/users/ash/BIN/genomics_general/VCF_processing/parseVCF.py \
-i Putative_Mutations.vcf.gz > Putative_Mutations.genotypes.txt
python /proj/snic2020-2-19/private/shark/users/ash/BIN/genomics_general/VCF_processing/parseVCF.py \
-i Putative_Mutations.AB-Filtered.vcf.gz > Putative_Mutations.AB-Filtered.genotypes.txt

#Extract allele counts for each genotype from VCF
bcftools query -f '%CHROM %POS[\t%AD]\n' Putative_Mutations.vcf.gz \
> Putative_Mutations.AD.txt

#For each putative mutation extract low quality genotypes from full dataset "geno" file created earlier
for MUTATION in {1..24}
do
POS=$(head -n $MUTATION mutation.list | tail -n 1)
cat Merged.sorted.noINDELS.HighlyMappable.NonRepeat.biallelic.genotypes.txt \
| grep "$POS" >> Putative_Mutations.missing.genotypes.txt
done 





##########################################################################
##########################################################################

#!/bin/bash
#SBATCH -A snic2021-5-8
#SBATCH -p core -N 1
#SBATCH -t 1-00:00:00
#SBATCH -J mpileupCalls

#Activate conda environment and load required modules
conda activate GenomicsGeneral
module load bioinfo-tools bcftools samtools/1.5

#For each putative mutation reconduct SNP calling using bcftools mpileup
REF=/proj/snic2020-2-19/private/shark/reference/satsuma/sHemOce1.mat.decon.20210528.fasta
mkdir bams
mkdir mpileup_VCFs

#First extract mutation positions from VCF
bcftools query -f '%CHROM\t%POS\n' ../identify_putative_mutations_gatk/Putative_Mutations.vcf.gz \
> mutation.positions

#Get list of samples
bcftools query -l ../identify_putative_mutations_gatk/Putative_Mutations.vcf.gz \
> sample.list

for MUTATION_NO in $(seq 1 $(wc -l < mutation.positions))
do
    SCAFF=$(head -n $MUTATION_NO mutation.positions | tail -n 1 | cut -f 1)
    POS=$(head -n $MUTATION_NO mutation.positions | tail -n 1 | cut -f 2)
    UPSTREAM=$(expr $POS - 5000)
    DOWNSTREAM=$(expr $POS + 5000)

    for SAMPLE in $(cat sample.list)
    do
        BAMS=/proj/snic2020-2-19/private/shark/users/mats/genotyping/${SAMPLE}_bam.list

        #Count number of BAM files in list
        BAM_COUNT=$(wc -l < $BAMS)

        #For each bam extract reads within region of interest
        for BAM_NO in $(seq 1 $BAM_COUNT)
        do
            samtools view -h $(head -n $BAM_NO $BAMS | tail -n 1 | cut -d " " -f 2) \
            "${SCAFF}:${UPSTREAM}-${DOWNSTREAM}" > temp_${BAM_NO}.bam
        done

        #Merge temp bam files
        ls temp_*.bam > temp.bam.list
        samtools merge -b temp.bam.list merged.bam -p

        #Fix sample names in header
        samtools view -H merged.bam | sed "s/SM:[^\t]*/SM:$SAMPLE/g" \
        | samtools reheader - merged.bam \
        > bams/${SCAFF}_${POS}_${SAMPLE}.bam

        #Tidy up
        rm merged.bam temp*
    done
    
    #Create list of sample bams for region of interest
    ls bams/${SCAFF}_${POS}_*.bam > bam.list

    #Re-call snps within region of interest using bcftools mpileup
    bcftools mpileup -Ou -f $REF -b bam.list -q 58 -Q 20 | bcftools call -m -Ov -o mpileup_VCFs/${SCAFF}_${POS}.vcf

    #Tidy up
    rm bam.list
done


##########################################################################
##########################################################################

#Pull out mutations from VCFs and comibine into a single new vcf file
#Extract header lines
cat $(ls mpileup_VCFs/scaffold_* | head -n 1) | grep "#" > putative_mutations.mpileup.vcf

#Extract mutation mpileup genotypes if present
for MUTATION_NO in $(seq 1 $(wc -l < mutation.positions))
do  
    SCAFF=$(head -n $MUTATION_NO mutation.positions | tail -n 1 | cut -f 1)
    POS=$(head -n $MUTATION_NO mutation.positions | tail -n 1 | cut -f 2)
    cat mpileup_VCFs/${SCAFF}_${POS}.vcf | grep -v "#" | grep ${POS} >> putative_mutations.mpileup.vcf
    echo "Processing complete for mutation no." $MUTATION_NO
done

#Convert to geno format
python /proj/snic2020-2-19/private/shark/users/ash/BIN/genomics_general/VCF_processing/parseVCF.py \
-i putative_mutations.mpileup.vcf > putative_mutations.mpileup.genotypes.txt
