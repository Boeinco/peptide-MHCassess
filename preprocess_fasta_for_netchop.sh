sed 's/>.*/|/g' $1 > protein_sequence_with_pipe.fasta
perl -pe 's/\s+//g' protein_sequence_with_pipe.fasta > $2
