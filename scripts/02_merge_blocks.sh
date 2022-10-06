#!/bin/bash
#SBATCH -A snic2022-5-242
#SBATCH -p core -N 1
#SBATCH -t 2-00:00:00
#SBATCH -J merge_blocks

#Load required modules
module load bioinfo-tools GATK/4.2.0.0 bcftools

#Merge blocks into a single VCF
cd joint_genotyping
ls block_*[0-9].vcf.gz > block.vcf.list
bcftools concat -f block.vcf.list -o raw_merged.vcf.gz -O z
gatk IndexFeatureFile -I raw_merged.vcf.gz
rm block_*[0-9].vcf.gz*
