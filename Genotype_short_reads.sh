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
