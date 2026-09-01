#!/bin/sh
#SBATCH --cpus-per-task=20
#SBATCH --job-name=eur_map
#SBATCH --partition=normal
#SBATCH --mem-per-cpu=10G


gatk GenomicsDBImport --java-options "-Xmx200G" \
    -L intervals.list \
    -R /biodata/dep_mercier/grp_novikova/A.Lyrata/ref/Alyrata/NT1_assembly/final_NT1/NT1_220222.fasta \
    -V samples.list \
 --overwrite-existing-genomicsdb-workspace \
 --reader-threads 25  --genomicsdb-workspace-path Pollen_AL62NT1_2n
 
 gatk --java-options "-Xmx85G -Xms78G -XX:ParallelGCThreads=8" GenotypeGVCFs -R NT1_220222.fasta -V gendb://Pollen_NT1NT8 --include-non-variant-sites -O Pollen_NT1NT8.vcf.gz

vcf="Pollen_NT1NT8_2n.vcf.gz"


# GATK best practices quality filter
bcftools filter -Ou -S . -e 'FMT/DP<4' ${vcf} |
  bcftools view -e "(N_MISSING==N_SAMPLES) || (QD<2.0 || FS>60.0 || SOR>3 || MQ<40.0 || MQRankSum<-12.0 || ReadPosRankSum<-8.0)" \
  --threads 128 -Oz -Wtbi -o ${vcf%.vcf*}_qual.vcf.gz

# Biallelic segregating sites
bcftools view -m2 -M2 -v snps -e 'AC==0 || AC==AN' ${vcf%.vcf*}_qual.vcf.gz -Oz -Wtbi --threads 128 -o ${vcf%.vcf*}_bial.vcf.gz
