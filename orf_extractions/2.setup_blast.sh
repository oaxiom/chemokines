gunzip -c orfs_BraLan.fa.gz   | makeblastdb -in - -out orfs_BraLan.db -dbtype prot -title orfs_BraLan
gunzip -c orfs_Eburgeri.fa.gz | makeblastdb -in - -out orfs_Eburgeri.db -dbtype prot -title orfs_Eburgeri
gunzip -c orfs_hg38.fa.gz     | makeblastdb -in - -out orfs_hg38.db -dbtype prot -title orfs_hg38
gunzip -c orfs_Pmarinus.fa.gz | makeblastdb -in - -out orfs_Pmarinus.db -dbtype prot -title orfs_Pmarinus