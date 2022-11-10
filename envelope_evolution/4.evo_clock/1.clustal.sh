cat ../3.envs_versus_orf_frags/envs.*.fasta >all_envs.fasta

clustalo --force -t Protein --threads=3 --outfmt=phylip -i all_envs.fasta -o all_envs.phylip  --guidetree-out=all_envs.dnd --wrap=120
