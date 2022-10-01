#!/bin/bash
#SBATCH -A snic2022-5-242
#SBATCH -p core -N 1
#SBATCH -t 4-00:00:00
#SBATCH --array=1-9:1
#SBATCH -J genotype_filtering

#Activate conda environment and load required modules
conda activate GenomicsGeneral
module load bioinfo-tools bcftools BEDTools/2.29.2 samtools/1.5 vcftools/0.1.16

#Set offspring ID using slurm_array_task_id
#OFFSPRING_ID=$(head -n $SLURM_ARRAY_TASK_ID ../resources/offspringIDs.txt | tail -n 1)
OFFSPRING_ID=ind1925

#Create / move into directory "genotype_filtering"
if [[ ! -d genotype_filtering ]]
then
    mkdir genotype_filtering
    mkdir genotype_filtering/${OFFSPRING_ID}
    cd genotype_filtering/${OFFSPRING_ID}
else
    mkdir genotype_filtering/${OFFSPRING_ID}
    cd genotype_filtering/${OFFSPRING_ID}
fi

#Called by GATK:
# - monomorphic / biallelic sites only
# - no indels
# - genotyped in offspring and both parents
#zcat ../joint_genotyping/trio_${OFFSPRING_ID}.sorted.vcf.gz \
zcat ../../joint_genotyping/block_1.trio_${OFFSPRING_ID}.vcf.gz \
| bcftools view --exclude-types indels --max-alleles 2 \
| bcftools view -e 'GT[0]="mis" || GT[1]="mis" || GT[2]="mis"' \
| grep -v "*" \
| bgzip > trio_${OFFSPRING_ID}.called_by_GATK.vcf.gz
COUNT=$(zcat trio_${OFFSPRING_ID}.called_by_GATK.vcf.gz | grep -v "#" | wc -l)
echo "Called by GATK: " $COUNT > trio_${OFFSPRING_ID}.site.counts.txt

#Remove repeat regions / regions of low mappability
bedtools intersect -a trio_${OFFSPRING_ID}.called_by_GATK.vcf.gz \
-b ../../resources/HiglyMappable_Unmasked_ranges.txt -header \
| bgzip > trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.vcf.gz
zcat trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.vcf.gz | grep -v "#" | wc -l \
> trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.counts.txt
COUNT=$(zcat trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.vcf.gz | grep -v "#" | wc -l)
echo "Non-repeat / highly mappable: " $COUNT >> trio_${OFFSPRING_ID}.site.counts.txt

#Remove genotypes with GQ < 20
bcftools query -l trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.vcf.gz \
> all.samples
zcat trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.vcf.gz \
| bcftools filter -S . -e 'FMT/GQ[@all.samples] < 20' \
| bcftools view -e 'GT[0]="mis" || GT[1]="mis" || GT[2]="mis"' \
| bgzip > trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.minGQ20.vcf.gz
COUNT=$(zcat trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.minGQ20.vcf.gz | grep -v "#" | wc -l)
echo "Min GQ20: " $COUNT >> trio_${OFFSPRING_ID}.site.counts.txt

###
#FILE STILL CONTAINS ASTERISKS IN ALT COLUMN -- THESE ARE DELETIONS (NEED TO REMOVE)

#Create a subset of snps where parents are homozygous for different alleles and all
#offspring are heterozygous and all sample genotypes are non-missing. These "known" heterozygous sites will
#be used as the basis for defining custom SNP filtering criteria.
zcat trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.minGQ20.vcf.gz \
| bcftools view --types snps \
| bcftools view -i 'GT[1]="ref" & GT[1]="hom"' \
| bcftools view -i 'GT[2]="alt" & GT[2]="hom"' \
| bcftools view -i 'GT[0]="het"' \
| bgzip > temp1.vcf.gz
zcat trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.minGQ20.vcf.gz \
| bcftools view --types snps \
| bcftools view -i 'GT[2]="ref" & GT[2]="hom"' \
| bcftools view -i 'GT[1]="alt" & GT[1]="hom"' \
| bcftools view -i 'GT[0]="het"' \
| bgzip > temp2.vcf.gz
bcftools index temp1.vcf.gz
bcftools index temp2.vcf.gz
bcftools concat --allow-overlaps temp1.vcf.gz temp2.vcf.gz \
| bgzip > trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.minGQ20.offspringHet.vcf.gz
rm temp1.vcf.gz*
rm temp2.vcf.gz*
COUNT=$(zcat trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.minGQ20.offspringHet.vcf.gz | grep -v "#" | wc -l)
echo "High confidence heterozygous sites: " $COUNT >> trio_${OFFSPRING_ID}.site.counts.txt

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
Rscript ../../scripts/generate.custom.filters.R \
trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.minGQ20.offspringHet.vcf.gz \
trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.minGQ20.vcf.gz \
trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.minGQ20.customFiltered.vcf.gz \
3

