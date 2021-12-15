makeblastdb -in simple_filtered.fasta -dbtype prot -title envs
blastp -query simple_filtered.fasta -db simple_filtered.fasta -evalue 1e-20 -num_threads 3 -outfmt "10 qaccver saccver pident length mismatch qstart qend sstart send evalue"  | uniq | gzip >simple_filtered.results.gz
