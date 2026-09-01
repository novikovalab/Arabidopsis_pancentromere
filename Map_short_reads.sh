#! /bin/bash
#SBATCH --array=1-4
#SBATCH --cpus-per-task=10
#SBATCH --job-name=eur_map
#SBATCH --partition=normal
#SBATCH --mem-per-cpu=5G
#SBATCH --output=job_output_%j.out
#SBATCH --error=job_output_%j.err
samples='down'
names='down'
###125

echo "STARTING JOB"
num=$SLURM_ARRAY_TASK_ID\p
acc=`sed -n ${num} ${names} | cut -d " " -f1`
name=`sed -n ${num} ${samples} | cut -d " " -f1`
#echo ${acc}
#acc="NT8_4_24"
#name="NT8_4_24"
#raw="/biodata/dep_mercier/grp_novikova/A.Lyrata/Robot25_lyrata_raw_novogene"
raw="."
out="."
####TRIMM READS

first=`ls ${name}*_1.fastq.gz`
second=`ls ${name}*_2.fastq.gz`
java -ea -Xmx18g -Xms1g -XX:ParallelGCThreads=4 -cp /opt/share/software/packages/BBMap_38.90/bin/current/ jgi.BBDuk t=4 in=${first} in2=${second} \
    out=${out}/${acc}_1.paired.fq.gz out2=${out}/${acc}_2.paired.fq.gz \
    outm=${out}/${acc}_1.unpaired.fq.gz outm2=${out}/${acc}_2.unpaired.fq.gz \
    ref=/netscratch/dep_mercier/grp_novikova/shared/adapters.fa \
    ktrim=r k=23 mink=11 hdist=1 tbo tpe qtrim=rl trimq=15 minlen=70
##############

ref='/biodata/dep_mercier/grp_novikova/A.Lyrata/ref/Alyrata/NT1_assembly/final_NT1/NT1_220222'
#raw="/biodata/dep_mercier/grp_novikova/A.Lyrata/lyrata_raw_data_November2021"
####bwa index ${ref}.fasta
bwa mem -t 4 -M -R '@RG\tID:lyr_'${acc}'\tSM:'${acc}'\tPL:Illumina\tLB:lyr_'${acc} ${ref}.fasta ${acc}_1.paired.fq.gz ${acc}_2.paired.fq.gz > ${out}/${acc}.sam
wait
#rm ${out}/${acc}_*.paired.fq
####samtools view -bh -F 0x2 ${out}/${acc}.sam > ${out}/${acc}.splitreads_filtered.bam
samtools view -bh -t ${ref}.fasta.fai -o ${out}/${acc}.bam ${out}/${acc}.sam
wait
#rm ${out}/${acc}.sam
samtools sort ${out}/${acc}.bam -o ${out}/${acc}.sort.bam -T ${out}/${acc}_temp
wait
rm ${out}/${acc}.bam
####samtools index ${out}/${acc}.splitreads_filtered.sort.bam
samtools sort -n ${out}/${acc}.sort.bam -o ${out}/${acc}.sortn.bam -T ${out}/${acc}_temp
wait
rm ${out}/${acc}.sort.bam
samtools fixmate -rm ${out}/${acc}.sortn.bam ${out}/${acc}.fixmate.sortn.bam
wait
rm ${out}/${acc}.sortn.bam
samtools sort ${out}/${acc}.fixmate.sortn.bam -o ${out}/${acc}.fixmate.sort.bam -T ${out}/${acc}_temp
wait
rm ${out}/${acc}.fixmate.sortn.bam
samtools markdup -rs ${out}/${acc}.fixmate.sort.bam ${out}/${acc}.fixmate.sort.markdup.bam
samtools index ${out}/${acc}.fixmate.sort.markdup.bam
wait
rm ${out}/${acc}.fixmate.sort.bam

gatk --java-options "-Xmx35G" HaplotypeCaller  \
   -R ${ref}.fasta \
   -I ${out}/${acc}.fixmate.sort.markdup.bam \
   -O ${out}/${acc}.g.vcf.gz \
   -ERC GVCF \
   --sample-ploidy 2 \
   --output-mode EMIT_ALL_CONFIDENT_SITES
