#!/bin/bash
#SBATCH -A snic2022-5-242
#SBATCH -p core -N 1
#SBATCH -t 2-00:00:00
#SBATCH --array=1-7:1
#SBATCH -J genotype_filtering

#Activate conda environment and load required modules
conda activate genomics_general 
module load bioinfo-tools GATK/4.2.0.0 bcftools BEDTools/2.29.2 samtools/1.5 vcftools/0.1.16

#Set offspring ID using slurm_array_task_id
OFFSPRING_ID=$(head -n $SLURM_ARRAY_TASK_ID resources/offspringIDs.txt | tail -n 1)

#Merge multiple block VCFs into a single VCF file
cd joint_genotyping
#ls block_*[0-9].trio_${OFFSPRING_ID}.vcf.gz > ${OFFSPRING_ID}.vcf.list
#bcftools concat -f ${OFFSPRING_ID}.vcf.list -o trio_${OFFSPRING_ID}.vcf.gz -O z
gatk IndexFeatureFile -I trio_${OFFSPRING_ID}.vcf.gz
rm block_*[0-9].trio_${OFFSPRING_ID}.vcf.gz

#Create / move into directory "genotype_filtering"
if [[ ! -d ../genotype_filtering ]]
then
    mkdir ../genotype_filtering
    mkdir ../genotype_filtering/${OFFSPRING_ID}
    cd ../genotype_filtering/${OFFSPRING_ID}
else
    mkdir ../genotype_filtering/${OFFSPRING_ID}
    cd ../genotype_filtering/${OFFSPRING_ID}
fi

#Called by GATK:
# - monomorphic / biallelic sites only
# - genotyped in offspring and both parents
RAW_VCF=../../joint_genotyping/trio_ind2024.vcf.gz
gatk SelectVariants \
     -R /proj/snic2020-2-19/private/shark/reference/satsuma/sHemOce1.mat.decon.20210528.fasta \
     -V $RAW_VCF \
     --select-type-to-include SNP \
     --select-type-to-include NO_VARIATION \
     -O temp.vcf.gz
bcftools view temp.vcf.gz --max-alleles 2 \
    | bcftools view -e 'GT[0]="mis" || GT[1]="mis" || GT[2]="mis"' \
    | bgzip > trio_${OFFSPRING_ID}.called_by_GATK.vcf.gz
rm temp.vcf.gz*
COUNT=$(zcat trio_${OFFSPRING_ID}.called_by_GATK.vcf.gz | grep -v "#" | wc -l)
echo "Biallelic / monomorphic called by GATK: " $COUNT > trio_${OFFSPRING_ID}.site.counts.txt

#Remove repeat regions / regions of low mappability
bedtools intersect -a trio_${OFFSPRING_ID}.called_by_GATK.vcf.gz \
-b ../../resources/HiglyMappable_Unmasked_ranges.txt -header \
| bgzip > trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.vcf.gz
COUNT=$(zcat trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.vcf.gz | grep -v "#" | wc -l)
echo "Non-repeat / highly mappable: " $COUNT >> trio_${OFFSPRING_ID}.site.counts.txt
rm trio_${OFFSPRING_ID}.called_by_GATK.vcf.gz*

#Remove genotypes with GQ < 20
bcftools query -l trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.vcf.gz \
> all.samples
zcat trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.vcf.gz \
| bcftools filter -S . -e 'FMT/GQ[@all.samples] < 20' \
| bcftools view -e 'GT[0]="mis" || GT[1]="mis" || GT[2]="mis"' \
| bgzip > trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.minGQ20.vcf.gz
COUNT=$(zcat trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.minGQ20.vcf.gz | grep -v "#" | wc -l)
echo "Min GQ20: " $COUNT >> trio_${OFFSPRING_ID}.site.counts.txt
rm trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.vcf.gz

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

#Based on the distribution of SNP quality annotations in the above VCF, define
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

#Apply filtering criteria
source apply.custom.filters.sh
COUNT=$(zcat trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.minGQ20.customFiltered.vcf.gz | grep -v "#" | wc -l)
echo "Callable sites: " $COUNT >> trio_${OFFSPRING_ID}.site.counts.txt

#From the custom filtered dataset, extract positons where both parents are homozygous
#for the REF allele and offsping is heterozygous
zcat trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.minGQ20.customFiltered.vcf.gz \
| bcftools view -i 'GT[1]="ref" & GT[1]="hom"' \
| bcftools view -i 'GT[2]="ref" & GT[2]="hom"' \
| bcftools view -i 'GT[0]="het"' \
| bgzip > temp1.vcf.gz

#From the custom filtered dataset, extract positons where both parents are homozygous
#for the ALT allele and offsping is heterozygous
zcat trio_${OFFSPRING_ID}.called_by_GATK.HighlyMappable.NonRepeat.minGQ20.customFiltered.vcf.gz \
| bcftools view -i 'GT[1]="alt" & GT[1]="hom"' \
| bcftools view -i 'GT[2]="alt" & GT[2]="hom"' \
| bcftools view -i 'GT[0]="het"' \
| bgzip > temp2.vcf.gz

#Combine the last two VCFs into a single file, these represent high confidence putative mutations
bcftools index temp1.vcf.gz
bcftools index temp2.vcf.gz
bcftools concat --allow-overlaps temp1.vcf.gz temp2.vcf.gz \
| bgzip > trio_${OFFSPRING_ID}.putative_mutations.vcf.gz
rm temp1.vcf.gz
rm temp2.vcf.gz
COUNT=$(zcat trio_${OFFSPRING_ID}.putative_mutations.vcf.gz | grep -v "#" | wc -l)
echo "Candidate DNMs: " $COUNT >> trio_${OFFSPRING_ID}.site.counts.txt

#Convert putative mutations to "geno" format
python /proj/snic2020-2-19/private/shark/users/ash/BIN/genomics_general/VCF_processing/parseVCF.py \
-i trio_${OFFSPRING_ID}.putative_mutations.vcf.gz > trio_${OFFSPRING_ID}.putative_mutations.genotypes.txt


###END