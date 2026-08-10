#This script chooses only bib scaffolds based on names. You can also choose only first n scaffolds.
from Bio import SeqIO

out=open('Pref.fa', 'w')
for seq_record in SeqIO.parse("Pref_rt/ragtag.scaffolds.fasta", 'fasta'):
    #["scaffold_1_RagTag", "scaffold_2_RagTag", "scaffold_3_RagTag", "scaffold_4_RagTag", "scaffold_5_RagTag", "scaffold_6_RagTag", "scaffold_7_RagTag", "scaffold_8_RagTag"]:
    if seq_record.id in ["scaffold_2_RagTag", "scaffold_1_RagTag", "scaffold_3_RagTag", "scaffold_4_RagTag", "scaffold_5_RagTag", "scaffold_6_RagTag", "scaffold_7_RagTag", "scaffold_8_RagTag"]:
    	SeqIO.write(seq_record, out, "fasta")
    #break

out.close()
