#PBS -l nodes=1:ppn=32
#PBS -j oe
#PBS -o results.txt
#PBS -q batch
#PBS -V 
cd $PBS_O_WORKDIR

blastp -max_target_seqs 10000 -query 1.blastdb/all_envs.fasta -db 1.blastdb/all_envs.fasta -evalue 1e-50 -num_threads 30 -outfmt "10 qaccver saccver pident length mismatch qstart qend sstart send evalue"  | gzip >all_vs_all.results.gz
