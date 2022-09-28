#!/bin/bash
#SBATCH -A snic2022-5-242
#SBATCH -p core -N 1
#SBATCH -t 4-00:00:00
#SBATCH --array=2-9:1
#SBATCH -J joint_genotyping

#Load required modules
module load bioinfo-tools GATK/4.2.0.0 bcftools

#Create / move into directory "joint_genotyping"
if [[ ! -d joint_genotyping ]]
then
    mkdir joint_genotyping
    cd joint_genotyping
else
    cd joint_genotyping
fi

#Set offspring ID using slurm_array_task_id
OFFSPRING_ID=$(head -n $SLURM_ARRAY_TASK_ID ../resources/offspringIDs.txt | tail -n 1)

#For each genome block do the following ...
for BLOCK_NO in $(seq 1 10)
do 
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

    #Remove gVCF
    rm block_${BLOCK_NO}.trio_${OFFSPRING_ID}.g.vcf.gz
done

#Merge multiple block VCFs into a single VCF file
#Create list of input VCFs
ls block_*[0-9].trio_${OFFSPRING_ID}.vcf.gz > ${OFFSPRING_ID}.vcf.list
bcftools concat -f ${OFFSPRING_ID}.vcf.list -o trio_${OFFSPRING_ID}.vcf.gz -O z
bcftools sort -o trio_${OFFSPRING_ID}.sorted.vcf.gz -O z trio_${OFFSPRING_ID}.vcf.gz
bcftools index trio_${OFFSPRING_ID}.sorted.vcf.gz

#Remove intermediate files
rm ${OFFSPRING_ID}.vcf.list
rm block_*[0-9].trio_${OFFSPRING_ID}.vcf.gz
rm trio_${OFFSPRING_ID}.vcf.gz
