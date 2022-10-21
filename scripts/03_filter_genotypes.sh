#!/bin/bash
#SBATCH -A snic2022-5-242
#SBATCH -p core -N 1
#SBATCH -t 1-00:00:00
#SBATCH -J FilterGenos

#Activate conda environment and load required modules
conda activate genomics_general
module load bioinfo-tools GATK/4.2.0.0 bcftools BEDTools/2.29.2 samtools/1.5 vcftools/0.1.16

#Make directory 
#mkdir genotype_filtering_redo
#cd genotype_filtering_redo

#From raw VCF file extract invariant sites
#plus remove any site that doesnt pass filtering (i.e. contains LowQual in filter column)
gatk SelectVariants \
    -R /proj/snic2020-2-19/private/shark/reference/satsuma/sHemOce1.mat.decon.20210528.fasta \
    -V ../joint_genotyping/raw_merged.vcf.gz \
    --select-type-to-include NO_VARIATION \
    -O called_by_GATK_invariant.vcf.gz

#From raw VCF file extract biallelic sites
#plus remove any site that doesnt pass filtering (i.e. contains LowQual in filter column)
gatk SelectVariants \
     -R /proj/snic2020-2-19/private/shark/reference/satsuma/sHemOce1.mat.decon.20210528.fasta \
     -V ../joint_genotyping/raw_merged.vcf.gz \
     --select-type-to-include SNP \
     --restrict-alleles-to BIALLELIC \
     -O called_by_GATK_biallelic.vcf.gz

#Combine invariant and biallelic VCFs
java -jar /proj/snic2020-2-19/private/shark/users/ash/BIN/picard/build/libs/picard.jar MergeVcfs \
    I=called_by_GATK_invariant.vcf.gz \
    I=called_by_GATK_biallelic.vcf.gz \
    O=called_by_GATK_invariant_plus_biallelic.vcf.gz

#Remove non-combined VCFs
rm called_by_GATK_invariant.vcf.gz rm called_by_GATK_biallelic.vcf.gz

#Remove sites where parental genotypes are missing as these aren't informative.
zcat called_by_GATK_invariant_plus_biallelic.vcf.gz \
| bcftools view -e 'GT[9]="mis" || GT[10]="mis"' \
| bgzip > called_by_GATK_invariant_plus_biallelic.informative_sites.vcf.gz

#Remove sites from repeat regions / regions of low mappability
bedtools intersect -a called_by_GATK_invariant_plus_biallelic.informative_sites.vcf.gz \
    -b ../resources/HiglyMappable_Unmasked_ranges.txt -header \
    | bgzip > called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.vcf.gz

#Set genotypes with GQ < 20 to missing "./."
bcftools query -l called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.vcf.gz > all.samples
zcat called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.vcf.gz \
| bcftools filter -S . -e 'FMT/GQ[@all.samples] < 20' \
| bcftools view -e 'GT[9]="mis" || GT[10]="mis"' \
| bgzip > called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.vcf.gz

#Set sample genotype to missing if allele balance is:
#< 0.2 for heterozygote calls
#< 0.9 for homozygote calls 
#This is implemented using the jvarkit tool "vcffilterjdk"
source ../scripts/filter_het_homo_by_AB.sh \
called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.vcf.gz \
temp.vcf.gz
zcat temp.vcf.gz | bcftools view -e 'GT[9]="mis" || GT[10]="mis"' \
| bgzip > called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.AB_filtered.vcf.gz
rm temp.vcf.gz

#Extract positions where parents are homozygous for different alleles
#and all offspring are heterozygous and all sample genotypes are non-missing.
zcat called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.AB_filtered.vcf.gz \
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
| bcftools sort | bgzip > temp1.vcf.gz

zcat called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.AB_filtered.vcf.gz \
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
| bcftools sort | bgzip > temp2.vcf.gz
bcftools index temp2.vcf.gz

java -jar /proj/snic2020-2-19/private/shark/users/ash/BIN/picard/build/libs/picard.jar MergeVcfs \
    I=temp1.vcf.gz \
    I=temp2.vcf.gz \
    O=HighConf_heterozygous_sites.vcf.gz

#Remove temp files
rm temp1.vcf.gz* temp2.vcf.gz*

#Based on the distribution of genotype quality annotations in the above VCF, define
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
Rscript ../scripts/generate.custom.filters.R \
HighConf_heterozygous_sites.vcf.gz \
called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.AB_filtered.vcf.gz \
called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.AB_filtered.PASSED_inhouse_filters.vcf.gz \
11

#Apply filtering criteria
source apply.custom.filters.sh

#From the custom filtered dataset, extract positons where at least one offspring is hereozygous
zcat called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.AB_filtered.PASSED_inhouse_filters.vcf.gz \
| bcftools view -i 'GT[9]="ref" & GT[9]="hom"' \
| bcftools view -i 'GT[10]="ref" & GT[10]="hom"' \
| bcftools view -i 'GT[0]="het" || GT[1]="het" || GT[2]="het" || GT[3]="het" || GT[4]="het" || GT[5]="het" || GT[6]="het" || GT[7]="het" || GT[8]="het"' \
| bcftools sort | bgzip > temp1.vcf.gz

#From the custom filtered dataset, extract positons where both parents are homozygous
#for the ALT allele and offsping is heterozygous
zcat called_by_GATK_invariant_plus_biallelic.informative_sites.HighlyMappable.NonRepeat.minGQ20.AB_filtered.PASSED_inhouse_filters.vcf.gz \
| bcftools view -i 'GT[9]="alt" & GT[9]="hom"' \
| bcftools view -i 'GT[10]="alt" & GT[10]="hom"' \
| bcftools view -i 'GT[0]="het" || GT[1]="het" || GT[2]="het" || GT[3]="het" || GT[4]="het" || GT[5]="het" || GT[6]="het" || GT[7]="het" || GT[8]="het"' \
| bcftools sort | bgzip > temp2.vcf.gz

#Combine the last two VCFs into a single file, these represent high confidence putative mutations
java -jar /proj/snic2020-2-19/private/shark/users/ash/BIN/picard/build/libs/picard.jar MergeVcfs \
    I=temp1.vcf.gz \
    I=temp2.vcf.gz \
    O=putative_mutations.vcf.gz

#Remove temp files
rm temp1.vcf.gz* temp2.vcf.gz*
