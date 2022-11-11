
"""

Extract the hits

"""

from glbase3 import *

def extract_orfs(blast_results_filename, orf_fasta_filename):
    species_name = blast_results_filename.replace('.txt', '').replace('results.', '')
    output_filename = f"envs.{species_name}.fasta"

    fastas_to_get = set([])

    blast_results = open(blast_results_filename, 'r')
    for lin in blast_results:
        lin = lin.strip().split(',')

        fastas_to_get.add(lin[1])
    blast_results.close()

    print(f'Found {len(fastas_to_get)} for {species_name}')

    fastas = utils.convertFASTAtoDict(f'../2.Genomes/{species_name}.fa.gz', gzip_input=True)
    # despite the name it really returns a list of dicts in a genelist-like structure.
    # rangle it to a proper dict:
    fastas = {i['name']: i['seq'].replace('-', '') for i in fastas}

    print(f'Loaded {len(fastas):,} ORFs')

    fasta_results = {}
    out = open(f'{output_filename}', 'w')
    for orf_name in fastas_to_get:
        out.write(f'>{orf_name}\n')
        out.write(f'{fastas[orf_name]}\n')
        fasta_results[orf_name] = fastas[orf_name]
    out.close()

    return fasta_results

if __name__ == '__main__':
    rp = extract_orfs('results.orfs_Pmarinus.txt', '')
    rh = extract_orfs('results.orfs_hg38.txt', '')
    re = extract_orfs('results.orfs_Eburgeri.txt', '')
    rb = extract_orfs('results.orfs_BraLan.txt', '')

    # Also save a merged FASTA
