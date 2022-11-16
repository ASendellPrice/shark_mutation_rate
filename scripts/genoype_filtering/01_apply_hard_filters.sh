#!/bin/bash
#SBATCH -A snic2022-5-242
#SBATCH -p core -N 1
#SBATCH -t 1-00:00:00
#SBATCH -J hard_filtering

######################################################################################
# This script was used to incrememntaly apply the following hard filters to GATk
# genotype calls:
# 1. Extract monomorphic & biallelic sites (exclude sites marked as low quality)
# 2. Remove non-informative sites (sites with missing parental genotypes)
# 3. Remove sites from repeat regions / regions with mappability < 1
# 4. Apply minumum genotype quality filter (GQ 20>)

# A. Sendell-Price, Nov 2022
######################################################################################

#Load required modules
module load bioinfo-tools GATK/4.2.0.0 BEDTools/2.29.2 bcftools 

#Make directory and move into it
mkdir hard_filtering
cd hard_filtering

#From raw VCF file extract invariant sites
gatk SelectVariants \
-V ../joint_genotyping/raw_merged.vcf.gz \
--select-type-to-include NO_VARIATION \
--exclude-filtered \
--sequence-dictionary /proj/snic2020-2-19/private/shark/reference/satsuma/sHemOce1.mat.decon.20210528.dict \
-O called_by_GATK_invariant.vcf.gz

#From raw VCF file extract biallelic sites
gatk SelectVariants \
-V ../joint_genotyping/raw_merged.vcf.gz \
--select-type-to-include SNP \
--restrict-alleles-to BIALLELIC \
--exclude-filtered \
--sequence-dictionary /proj/snic2020-2-19/private/shark/reference/satsuma/sHemOce1.mat.decon.20210528.dict \
-O called_by_GATK_biallelic.vcf.gz

#Combine invariant and biallelic VCFs
java -jar /proj/snic2020-2-19/private/shark/users/ash/BIN/picard/build/libs/picard.jar MergeVcfs \
I=called_by_GATK_invariant.vcf.gz \
I=called_by_GATK_biallelic.vcf.gz \
O=called_by_GATK_invariant_plus_biallelic.vcf.gz \
SEQUENCE_DICTIONARY=/proj/snic2020-2-19/private/shark/reference/satsuma/sHemOce1.mat.decon.20210528.dict

#Remove un-needed VCFs
#rm called_by_GATK_invariant.vcf.gz
#rm called_by_GATK_biallelic.vcf.gz

#Remove sites from repeat regions / regions of low mappability
bedtools intersect -a called_by_GATK_invariant_plus_biallelic.vcf.gz \
-b ../resources/HiglyMappable_Unmasked_ranges.txt -header \
| bgzip > called_by_GATK_invariant_plus_biallelic.HighlyMappable.NonRepeat.vcf.gz

#Set genotypes with GQ < 20 to missing "./."
bcftools query -l called_by_GATK_invariant_plus_biallelic.HighlyMappable.NonRepeat.vcf.gz > all.samples
zcat called_by_GATK_invariant_plus_biallelic.HighlyMappable.NonRepeat.vcf.gz \
| bcftools filter -S . -e 'FMT/GQ[@all.samples] < 20' \
| bcftools view -e 'GT[9]="mis" || GT[10]="mis"' \
| bgzip > called_by_GATK_invariant_plus_biallelic.HighlyMappable.NonRepeat.minGQ20.vcf.gz

#Remove sites where parental genotypes are missing as these aren't informative.
zcat called_by_GATK_invariant_plus_biallelic.HighlyMappable.NonRepeat.minGQ20.vcf.gz \
| bcftools view -e 'GT[9]="mis" || GT[10]="mis"' \
| bgzip > called_by_GATK_invariant_plus_biallelic.HighlyMappable.NonRepeat.minGQ20.parental_non_missing.vcf.gz
