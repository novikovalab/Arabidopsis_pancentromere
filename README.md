# Arabidopsis pancentromere
This repository contains scripts used for centromere analysis of Arabidopsis species.
# Assembly
```mabs-hifiasm.py --pacbio_hifi_reads Raw_reads.fastq.gz --output_folder Pref --local_busco_dataset brassicales_odb10 --threads 30 # used ploidy parameter for tetraploids
ragtag.py scaffold  -f 100000 NT1_v2.fa Pref/The_best_assembly/assembly.fasta -o Pref_rt
```
Then keep only chromosomes, not the small contigs using **Subset_only_chromosomes.py**.

## Run Helixer (gene annotation)
```helixer_v0.3.3 Helixer.py --fasta-path Pref.fa --gff-output-path Pref_helixer.gff --lineage land_plant
```
## GENESPACE
Use **Data_preparation_for_genespace.sh** then **run_genespace3.R**.

## Tidk
```tidk search --string TTTAGGG --output Pref --dir . Pref.fa
tidk plot -o Pref --tsv Pref_telomeric_repeat_windows.tsv
```

## TRASH
```
../trash/TRASH/TRASH_run.sh Pref.fa --seqt NT1_repeat_template.csv --par 20 --horclass 'pAge2' #horclass for A. arenosa - pAa
```

## EDTA
```EDTA.pl --genome Pref.fa --force 1 --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 40 --curatedlib TAIR10_TE.liban 
```

## To see if the assembly is good load Helixer, tidk, trash, EDTA, bam to IGV and check 
- small contigs with no good reads mapped (random repeats ragtaged between good contigs)
- especially often small contigs on the chromosome ends, maybe after telomere, no contigs after telomere should stay
- centromeres on the chromosome ends
- duplications obvious from GENESPACE or duplicated centromeres on different contigs.
### How to cut duplications?
1) remove the contig, if the whole thing is duplicated
2) if many duplications can try to run purge_dups, it was not too useful though. Can run several iterations.
3) remove part of the contig, if only beginning of the contig overlaps with the previous one – take first half of the chromosome (one duplicated copy), and 100kb from the beginning of the second half (second copy), make 2 bed files,
getfasta for both, minimap qu on ref, look at bam file – should have good (primary) hit from the beginning of query sequence. Cut here.

- other weird things – 2 chromosomes glued together, parts missing or bad assembled, misassembly with small duplication, inverted contig due to large inversion in it. 
- inversions of the whole contig, whole chromosome for collaborators’ genomes, we invert (in python), misassembly we cut by the point with no coverage.
Still bad things? Might want to look at hap1 and hap2 assemblies, sometimes one of the haplotypes assembled better. Reassemble with other version of hifiasm, choose the chromosome which looks the best similar as for the primary genome.

Now we identified all things we want to cut and put them in a bed file. Keep 100 Ns between contigs!

### How we cut things?

```
bedtools complement -i Pref_cut.bed -g Pref.fa.fai > Pref_keep.bed # If bed of what to cut out
bedtools getfasta -fi Pref.fa -bed Pref_keep.bed > Pref_cut.fa
```

Remove the chromosome names I don’t need with sed (like when you have 2 parts of the same chromosome one after another and you don’t need to flip them), remove coordinates and "RagTag".
```
sed -i 's/:.*//' Pref_cut.fa
```
Make multiline fasta and invert contigs if needed with SeqIO in python.
Then check the assembly again with GENESPACE, trash, tidk.

## Identify centromere boundaries

Open EDTA, trash, ragtag in IGV and look from the first to last centromere repeat. 
All trash centromeric repeats should be inside. EDTA annotates them as "repeat region 18". If they cover more than trash, cut by them unless other TE overlaps, then by this TE but always keep all trash inside.

## Extract repeat libraries from TRASH results
Filter centromeres for PCAs
```tail -n +2 all.repeats.from.Pref.fa.csv | awk -F "," '{print $8, $1, $2, $3, $4}' | sed 's/\"//g' | sed 's/ /\t/g'> all.repeats.from.Pref.bed
bedtools intersect -a all.repeats.from.Pref.bed -wa -b Pref_centromeres.bed > Pref_centromere_repeats.bed
cut -f 4 Pref_centromere_repeats.bed | sort | uniq -c | less # To identify frequent repeat sizes and filter out too long or too short. The boundaries individually for each genome. At least 100 repeats present of this repeat length.
#[138/367)
awk '{if ($4 > 138) print $1, $2, $3, $4, $5}' Pref_centromere_repeats.bed | awk '{if ($4 < 367) print $1, $2, $3, $4, $5}' > Pref_centromere_filt_repeats.bed
```
Then PCA script in python to classify the repeats (**PCA_script.py**).
Reclassify if not enough pAge2, or A. arenosa, or remove 4th repeat type and reclassify rest for A. cebennensis and A. pedemontana (4 repeat types are too much for one PCA, and the types are getting mixed).
The PCA reclassification used for the centromere repeat plot (**Big_centromere_repeat_plot.py**).

## Making trees
Subsample repeats
```bedtools getfasta -bed <(grep pAa Repeats_with_classes.bed | shuf -n 10 ) -fi Pref.fa | sed 's/>/\>Pref_pAa/' >> Subset_repeats.fa
#Align:
linsi --adjustdirection Subset_repeats.fa > Subset_repeats_align.fa
```
Here open in Ugene and trim tails that are not aligned properly and remove any sequences which make no sense (not aligned well).
Tree:
```iqtree -bb 1000 -s Split_pAa_and_CEN178_fin1.fa
```
Vizualize in Itol-tree.

## BLAST
Blasting of both canonical and non-canonical repeats to know where they are on chromosome arms (fig. S13, S15).
```blastn -db Pref.fa -query Reference_rep.fa -outfmt "6 sseqid sstart send qseqid bitscore sseq" >> Pref_blast_rep.bed
bedtools intersect -v -a <(awk '{if ($2 < $3) print $0; else print $1,$3,$2,$4,$5}' Pref_blast_rep.bed | sed 's/ /\t/g') -b Pref_centromeres.bed > Pref_canonical_nocentr.bed
#Second level of blast – blast something from outside centromeres to centromeres
#Get fastas like this: 
bedtools getfasta -bed <(cat Pref_canonical_nocentr.bed | awk '{if ($3 < $2) print $1, $3, $2; else print $1,$2,$3}' | sed 's/ /\t/g') -fi Pref.fa
while IFS=" " read -r ID remainder
do
blastn -db all_centromeres_only.fa -query ${ID}_nocentr.fa -outfmt "6 sseqid sstart send qseqid bitscore" -max_hsps 1 -subject_besthit -max_target_seqs 1 -num_threads 20 > ${ID}_can_blasted_only 
done < "listtoblast"
```
## 60-kmers analysis 
```for samp in `ls all_centromeres_only.fa.split`; do
jellyfish count -m 60 -s 100M all_centromeres_only.fa.split/$samp -o splitcen/$samp
jellyfish dump -c splitcen/$samp | sort -k2,2rn > splitcen/$samp.counts
done
```
For A. kamchatica k-mers (fig. S20):
chose 500 to 800 line when sorted by freqeuncy, then:
```for kmer in `cat jf_kron_chr5_pAge1_kmers`; do
seqkit locate --pattern $kmer  all_centromeres.fa | cut -f 1 | uniq  > kmer_stuff
wc -l kmer_stuff >> kron_chr5pAge1_hits
#filtering out something not specific, with hits on more than 51 scaffolds
if (( $(wc -l < kmer_stuff) < 51)); then
cat kmer_stuff >> all_kmers_chr5_pAge1_kron
fi
done
```
## Simulations
We used scripts from [https://github.com/schneebergerlab/replicated-assemblies-centromere-study]: https://github.com/schneebergerlab/replicated-assemblies-centromere-study
We also added blasting on every step to break the simulation if no inserted repeats left. A sample script **Centromere_simulation.sh** has dependances in replicated-assemblies-centromere-study repository. 

## Short reads map to centromere
The mapping performed with **map_short_reads_on_repeats.sh**, then map vizualization done in QGIS.

## Read mapping and genotype calling
Was used for QTL and for pollen sequencing data. The mapping script: **Map_short_reads.sh**, combining, genotyping and filtering: **Genotype_short_reads.sh**.


## Pollen sequencing analysis
Make ped files for every leaf and pollen F1 sample. Example ped file:
X04     SRR8157565      SRR8157563      SRR8157564      0       1

For leaf and pollen run this:
```#MIND THE ORDER
vcf-subset --exclude-ref -c 25114-2,24030-4,24036-1 Pollen_AL62NT1_2n_bial.vcf.gz > X02_L.vcf.gz &
whatshap phase --ped arabidopsis_SD_workflow/SD_phasing_workflow/workflow/peds/X01_L.ped --ref=/biodata/dep_mercier/grp_novikova/A.Lyrata/ref/Alyrata/NT1_assembly/final_NT1/NT1_220222.fasta --recombrate 4.8  -o X01_L.phased.vcf X01_L.vcf.gz 24030-4.fixmate.sort.markdup.bam 24036-1.fixmate.sort.markdup.bam /biodata/dep_mercier/grp_novikova/A.Lyrata/Robot25_lyrata_raw_novogene/bams/25114-1.fixmate.sort.markdup.bam &
python arabidopsis_SD_workflow/SD_phasing_workflow/workflow/scripts/filter_vcf_with_stats.py Control_L.phased.vcf Control_P.phased.vcf Control_L Control_P Control_L.filtered.vcf Control_P.filtered.vcf
```

Then we need to make a mpileup from only forward reads:
```python arabidopsis_SD_workflow/SD_phasing_workflow/workflow/scripts/trim_reads_v4.py SRR8157566_1.paired.fq.gz SRR8157565_1.paired.fq.gz SRR8157566 SRR8157565 ./ &
bwa mem -t 4 /biodata/dep_mercier/grp_novikova/A.Lyrata/ref/Alyrata/NT1_assembly/final_NT1/NT1_220222.fasta AL62TE11_pollen.retrimmed.fastq | samtools view -Sb | samtools sort > AL62TE11_1_pollen.retrimmed.bam &
bcftools mpileup -a AD --fasta-ref /biodata/dep_mercier/grp_novikova/A.Lyrata/ref/Alyrata/NT1_assembly/final_NT1/NT1_220222.fasta 25114-1.retrimmed.bam AL62TE11_1_pollen.retrimmed.bam > X01_2.mpileup &
bcftools call -mv -Ov X01_2.mpileup | bcftools norm -m -any | bcftools filter -e 'QUAL<20 || DP<2' > x01_mpileup/X01.vcf
mv x01_mpileup/X01.vcf x01_mpileup/X01_L.vcf
bgzip X01_L.vcf #and keep only one gz file in the folder
```
Now we can proceed.
try2.map is a decoy map looking like this:
scaffold_1      1000    0
scaffold_1      28694250        0
```# change index manually in the script, I also made some changes in the script
python arabidopsis_SD_workflow/SD_phasing_workflow/workflow/scripts/update_filt_vcf_depths.py X01_L.filtered.vcf X01_L x01_mpileup/

python ./scripts/extract_haps_v2.py X01_L.update.vcf X01_P.update.vcf  try2.map > X01.update.haplotypes
#need to split per chromosome otherwise the C++ script crashes
grep ^7 X02.final.haplotypes > X02.chr7.haplo
MAP_SD/src/try_map -d X02.chr7.haplo -w 10000 > X02.chr7.window_10K.likelihoods

# now put everything back again, or doesn’t work
cat X02.chr?.window_10K.likelihoods > X02.allmerge.window_10K.likelihoods
cat X02.chr?.haplo > X02.haplotypes
python workflow/scripts/plot_haps_v3.py X02.haplotypes X02.allmerge.window_10K.likelihoods X02.ratios.png
```
## SweeD
```SweeD  -name WestSib -outgroup Acebennensis1 -input All_lyrata_final_allpos_bial_dipl.vcf -osf Wsib.sf -sampleList ../../Eur_lyrata/Wsib_list
/netscratch/dep_mercier/grp_novikova/mvasilarou/bin/sweed/SweeD -name Wsib_grid500 -input Wsib.sf -grid 500
```

# Map and genotype long reads
Mapping with winnowmap:
```meryl count k=15 output merylDB_te8 Pref.fa
meryl print greater-than distinct=0.9998 merylDB_pref > pref_k15.txt
winnowmap -W pref_k15.txt -ax map-pb Pref.fa Reads.fastq.gz | samtools view -Sb | samtools sort > Pref_winnow_back.bam
```
The same way we mapped to NT1 reference genome. Unforunately, most software options don't work for mapping autopolyploids. So, we use GATK and need to add read groups.

```picard AddOrReplaceReadGroups  -I AL06_on_NT1.bam  -O AL08_to_NT1_rg.bam  --RGID 4  --RGLB lib1  --RGPL PACBIO  --RGPU unit1  --RGSM 20 --VALIDATION_STRINGENCY LENIENT
# then GATK with --includeNonVariantSites option
```
