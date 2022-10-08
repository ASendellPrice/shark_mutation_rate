#!/bin/bash
#SBATCH -A snic2022-5-242
#SBATCH -p core -N 1
#SBATCH -t 1-00:00:00
#SBATCH -J FilterGenos

#Activate conda environment and load required modules
conda activate GenomicsGeneral
module load bioinfo-tools GATK/4.2.0.0 bcftools BEDTools/2.29.2 samtools/1.5 vcftools/0.1.16

#Make directory 
mkdir genotype_filtering
cd genotype_filtering

#Called by GATK:
# - monomorphic / biallelic sites only
# - genotyped in offspring and both parents
gatk SelectVariants \
     -R /proj/snic2020-2-19/private/shark/reference/satsuma/sHemOce1.mat.decon.20210528.fasta \
     -V ../joint_genotyping/raw_merged.vcf.gz \
     --select-type-to-include SNP \
     --select-type-to-include NO_VARIATION \
     -O temp.vcf.gz
bcftools view temp.vcf.gz --max-alleles 2 \
    | bgzip > called_by_GATK.vcf.gz
rm temp.vcf.gz*

#Remove repeat regions / regions of low mappability
bedtools intersect -a called_by_GATK.vcf.gz \
-b ../resources/HiglyMappable_Unmasked_ranges.txt -header \
| bgzip > called_by_GATK.HighlyMappable.NonRepeat.vcf.gz

#Set genotypes with GQ < 20 to missing "./."
bcftools query -l called_by_GATK.HighlyMappable.NonRepeat.vcf.gz \
> all.samples
zcat called_by_GATK.HighlyMappable.NonRepeat.vcf.gz \
| bcftools filter -S . -e 'FMT/GQ[@all.samples] < 20' \
| bgzip > called_by_GATK.HighlyMappable.NonRepeat.minGQ20.vcf.gz

#Create a subset of snps where parents are homozygous for different alleles and all
#offspring are heterozygous and all sample genotypes are non-missing. These "known" heterozygous sites will
#be used as the basis for defining custom SNP filtering criteria.
zcat called_by_GATK.HighlyMappable.NonRepeat.minGQ20.vcf.gz \
| bcftools filter -e 'GT[@all.samples]="mis"' \
| bcftools filter -i 'GT[0]="het"' \
| bcftools filter -i 'GT[1]="het"' \
| bcftools filter -i 'GT[2]="het"' \
| bcftools filter -i 'GT[3]="het"' \
| bcftools filter -i 'GT[4]="het"' \
| bcftools filter -i 'GT[5]="het"' \
| bcftools filter -i 'GT[6]="het"' \
| bcftools filter -i 'GT[7]="het"' \
| bcftools filter -i 'GT[8]="het"' \
| bcftools view -i 'GT[9]="ref" & GT[9]="hom"' \
| bcftools view -i 'GT[10]="alt" & GT[10]="hom"' \
| bgzip > temp1.vcf.gz
bcftools index temp1.vcf.gz

zcat called_by_GATK.HighlyMappable.NonRepeat.minGQ20.vcf.gz \
| bcftools filter -e 'GT[@all.samples]="mis"' \
| bcftools filter -i 'GT[0]="het"' \
| bcftools filter -i 'GT[1]="het"' \
| bcftools filter -i 'GT[2]="het"' \
| bcftools filter -i 'GT[3]="het"' \
| bcftools filter -i 'GT[4]="het"' \
| bcftools filter -i 'GT[5]="het"' \
| bcftools filter -i 'GT[6]="het"' \
| bcftools filter -i 'GT[7]="het"' \
| bcftools filter -i 'GT[8]="het"' \
| bcftools view -i 'GT[9]="alt" & GT[9]="hom"' \
| bcftools view -i 'GT[10]="ref" & GT[10]="hom"' \
| bgzip > temp2.vcf.gz
bcftools index temp2.vcf.gz

bcftools concat --allow-overlaps temp1.vcf.gz temp2.vcf.gz \
| bgzip > called_by_GATK.HighlyMappable.NonRepeat.minGQ20.offspringHet.vcf.gz
rm temp1.vcf.gz*
rm temp2.vcf.gz*

#Based on the distribution of SNP quality annotations in the above VCF, define
#the filtering criteria (lower = 5th percentile, upper = 95th percentile).
#Filtering will be applied to the following site annotations:
#mapping quality: MQ (lower bound only)
#mapping quality rank sum: MQRankSum
#base quality rank sum: BaseQRankSum
#read position rank sum: ReadPosRankSum
#quality by depth: QD (lower bound only)
#Filtering will then be applied to the following individual annotations:
#genotype quality: GQ (lower bound only)
#sample read depth: DP 
Rscript ../../scripts/generate.custom.filters.R \
called_by_GATK.HighlyMappable.NonRepeat.minGQ20.offspringHet.vcf.gz \
called_by_GATK.HighlyMappable.NonRepeat.minGQ20.vcf.gz \
called_by_GATK.HighlyMappable.NonRepeat.minGQ20.customFiltered.vcf.gz \
11

#Apply filtering criteria
source apply.custom.filters.sh

#From the custom filtered dataset, extract positons where both parents are homozygous
#for the REF allele and at least one offspring is hereozygous
zcat called_by_GATK.HighlyMappable.NonRepeat.minGQ20.customFiltered.vcf.gz \
| bcftools view -i 'GT[9]="ref" & GT[9]="hom"' \
| bcftools view -i 'GT[10]="ref" & GT[10]="hom"' \
| bcftools view -i 'GT[0-8]="het"' \
| bgzip > temp1.vcf.gz

#From the custom filtered dataset, extract positons where both parents are homozygous
#for the ALT allele and offsping is heterozygous
zcat called_by_GATK.HighlyMappable.NonRepeat.minGQ20.customFiltered.vcf.gz \
| bcftools view -i 'GT[9]="alt" & GT[9]="hom"' \
| bcftools view -i 'GT[10]="alt" & GT[10]="hom"' \
| bcftools view -i 'GT[0-8]="het"' \
| bgzip > temp2.vcf.gz

#Combine the last two VCFs into a single file, these represent high confidence putative mutations
bcftools index temp1.vcf.gz
bcftools index temp2.vcf.gz
bcftools concat --allow-overlaps temp1.vcf.gz temp2.vcf.gz \
| bgzip > putative_mutations.vcf.gz
rm temp1.vcf.gz*
rm temp2.vcf.gz*
