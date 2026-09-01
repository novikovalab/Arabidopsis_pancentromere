ins="100"
for run_num in {41..50}; do
python ./bin/simulate.py \
    ./thaliana_ins${ins}.fa 1000 \
    pos_ins${ins}.pos \
    fasta/1000generation_ins${ins}_${run_num}.out.fa \
    records/1000generation_ins${ins}_${run_num}.record.txt \
    generation/1000generation_ins${ins}_${run_num}.unit.pos
echo ">ins1_1000gen" | cat - fasta/1000generation_ins${ins}_${run_num}.out.fa \
            > 1.fasta/1000generation_ins${ins}_${run_num}.out.fa
makeblastdb -in 1.fasta/1000generation_ins${ins}_${run_num}.out.fa -dbtype nucl
blastn -db 1.fasta/1000generation_ins${ins}_${run_num}.out.fa -query ../BaltAre_pAa_cons.fa -outfmt 6 | wc -l > sim_${ins}_${run_num}

for i in $(seq 2000 1000 1000000); do
    prev=$((i - 1000))
    python ./bin/simulate.py \
        fasta/${prev}generation_ins${ins}_${run_num}.out.fa 1000 \
        generation/${prev}generation_ins${ins}_${run_num}.unit.pos \
        fasta/${i}generation_ins${ins}_${run_num}.out.fa \
        records/${i}generation_ins${ins}_${run_num}.record.txt \
        generation/${i}generation_ins${ins}_${run_num}.unit.pos
        echo ">ins1_${i}gen" | cat - fasta/${i}generation_ins${ins}_${run_num}.out.fa \
            > 1.fasta/${i}generation_ins${ins}_${run_num}.out.fa
        makeblastdb -in 1.fasta/${i}generation_ins${ins}_${run_num}.out.fa -dbtype nucl
        blastn -db 1.fasta/${i}generation_ins${ins}_${run_num}.out.fa -query ../BaltAre_pAa_cons.fa -outfmt 6 | wc -l >> sim_${ins}_${run_num}
        num_paa=$(blastn -db 1.fasta/${i}generation_ins${ins}_${run_num}.out.fa -query ../BaltAre_pAa_cons.fa -outfmt 6 | wc -l)
        echo ${num_paa}
        if [[ "$num_paa" -eq 0 ]]; then
    break
fi
  done

done
