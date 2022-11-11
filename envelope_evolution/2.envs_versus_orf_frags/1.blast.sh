query='../1.nonredundant_envs/simple_filtered.fasta'
outfmt="\'10 qaccver saccver pident length mismatch qstart qend sstart send evalue\'"

# Can't put option in an opt for some reason.
blastp -query $query -db ../../orf_extractions/orfs_BraLan.db -evalue 1e-2 -num_threads 3  -outfmt '10 qaccver saccver pident length mismatch qstart qend sstart send evalue' >results.orfs_BraLan.txt
blastp -query $query -db ../../orf_extractions/orfs_Eburgeri.db -evalue 1e-2 -num_threads 3 -outfmt '10 qaccver saccver pident length mismatch qstart qend sstart send evalue' >results.orfs_Eburgeri.txt
blastp -query $query -db ../../orf_extractions/orfs_hg38.db -evalue 1e-2 -num_threads 3     -outfmt '10 qaccver saccver pident length mismatch qstart qend sstart send evalue' >results.orfs_hg38.txt
blastp -query $query -db ../../orf_extractions/orfs_Pmarinus.db -evalue 1e-2 -num_threads 3 -outfmt '10 qaccver saccver pident length mismatch qstart qend sstart send evalue' >results.orfs_Pmarinus.txt
