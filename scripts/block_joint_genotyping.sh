#!/bin/bash
#SBATCH -A snic2022-5-242
#SBATCH -p core -N 1
#SBATCH -t 2-00:00:00
#SBATCH --array=1-10:1
#SBATCH -J joint_genotyping

#Submitted like so:
#mkdir joint_genotyping
#cd joint_genotyping
#for OFFSPRING_ID in $(cat ../resources/offspringIDs.txt)
#do
#    sbatch ../scripts/block_joint_genotyping.sh $OFFSPRING_ID
#done

#Load required modules
module load bioinfo-tools GATK/4.2.0.0 bcftools

#Set offspring ID
OFFSPRING_ID=$1

#Set block number
BLOCK_NO=$SLURM_ARRAY_TASK_ID

#Combine GVCFs
gatk --java-options "-Xmx120g" CombineGVCFs \
-R /proj/snic2020-2-19/private/shark/reference/satsuma/sHemOce1.mat.decon.20210528.fasta \
--variant /proj/snic2020-2-19/private/shark/variants/gvcfs/${OFFSPRING_ID}_block_${BLOCK_NO}.g.vcf.gz \
--variant /proj/snic2020-2-19/private/shark/variants/gvcfs/sHemOce2_block_${BLOCK_NO}.g.vcf.gz \
--variant /proj/snic2020-2-19/private/shark/variants/gvcfs/sHemOce3_block_${BLOCK_NO}.g.vcf.gz \
--convert-to-base-pair-resolution \
-O block_${BLOCK_NO}.trio_${OFFSPRING_ID}.g.vcf.gz

#Conduct joint genotyping using gVCFs
gatk --java-options "-Xmx120g" GenotypeGVCFs \
--include-non-variant-sites \
-R /proj/snic2020-2-19/private/shark/reference/satsuma/sHemOce1.mat.decon.20210528.fasta \
-V block_${BLOCK_NO}.trio_${OFFSPRING_ID}.g.vcf.gz \
-O block_${BLOCK_NO}.trio_${OFFSPRING_ID}.vcf.gz

#Remove gVCF and its index
rm block_${BLOCK_NO}.trio_${OFFSPRING_ID}.g.vcf.gz*
