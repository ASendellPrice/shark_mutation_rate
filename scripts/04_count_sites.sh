#!/bin/bash
#SBATCH -A snic2022-5-242
#SBATCH -p core -N 1
#SBATCH -t 2-00:00:00
#SBATCH --array=1-9:1
#SBATCH -J count_sites

#load required modules
module load bioinfo-tools bcftools 

#Set offspring ID using slurm_array_task_id
OFFSPRING_ID=$(head -n $SLURM_ARRAY_TASK_ID resources/offspringIDs.txt | tail -n 1)

#Make directory for trio and move into it
mkdir genotype_filtering/${OFFSPRING_ID}_site_counts 
cd genotype_filtering/${OFFSPRING_ID}_site_counts 

#Make  a list of sample names for trio
echo $OFFSPRING_ID > trio.samples
tail -n 2 ../all.samples >> trio.samples

#Count sites in "called_by_GATK.vcf.gz"
# - monomorphic / biallelic sites only
# - genotyped in offspring and both parents
COUNT=$(bcftools view -e 'FILTER="LowQual"' ../called_by_GATK.vcf.gz | bcftools filter -e 'GT[@trio.samples]="mis"' | grep -v "#" | wc -l)
echo "Biallelic / monomorphic called by GATK: " $COUNT > trio_${OFFSPRING_ID}.site.counts.txt

#Count sites in "called_by_GATK.HighlyMappable.NonRepeat.vcf.gz"
COUNT=$(bcftools view -e 'FILTER="LowQual"' ../called_by_GATK.HighlyMappable.NonRepeat.vcf.gz | bcftools filter -e 'GT[@trio.samples]="mis"' | grep -v "#" | wc -l)
echo "Non-repeat / highly mappable: " $COUNT > trio_${OFFSPRING_ID}.site.counts.txt

#Count sites in "called_by_GATK.HighlyMappable.NonRepeat.minGQ20.vcf.gz"
COUNT=$(bcftools view -e 'FILTER="LowQual"' ../called_by_GATK.HighlyMappable.NonRepeat.minGQ20.vcf.gz | bcftools filter -e 'GT[@trio.samples]="mis"' | grep -v "#" | wc -l)
echo "Min GQ20: " $COUNT >> trio_${OFFSPRING_ID}.site.counts.txt

#Count sites in "called_by_GATK.HighlyMappable.NonRepeat.minGQ20.customFiltered.vcf.gz"
COUNT=$(bcftools view -e 'FILTER="LowQual"' ../called_by_GATK.HighlyMappable.NonRepeat.minGQ20.customFiltered.vcf.gz | bcftools filter -e 'GT[@trio.samples]="mis"' | grep -v "#" | wc -l)
echo "Callable sites: " $COUNT >> trio_${OFFSPRING_ID}.site.counts.txt
