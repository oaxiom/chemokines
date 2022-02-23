# extract all possible ORFs in all reading frames for some set of length criteria;

import sys, os, gzip, re

# Minimum length of ORF to consider
min_length_to_consider = 50

rcd = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

def rc(seq):
    return "".join(reversed([rcd[i] for i in seq]))

codons = {
    "AAA" : "K",    "AAC" : "N",    "AAG" : "K",    "AAT" : "N",
    "ACA" : "T",    "ACC" : "T",    "ACG" : "T",    "ACT" : "T",
    "AGA" : "R",    "AGC" : "S",    "AGG" : "R",    "AGT" : "S",
    "ATA" : "I",    "ATC" : "I",    "ATG" : "M",    "ATT" : "I",
    "CAA" : "Q",    "CAC" : "H",    "CAG" : "Q",    "CAT" : "H",
    "CCA" : "P",    "CCC" : "P",    "CCG" : "P",    "CCT" : "P",
    "CGA" : "R",    "CGC" : "R",    "CGG" : "R",    "CGT" : "R",
    "CTA" : "L",    "CTC" : "L",    "CTG" : "L",    "CTT" : "L",
    "GAA" : "E",    "GAC" : "D",    "GAG" : "E",    "GAT" : "D",
    "GCA" : "A",    "GCC" : "A",    "GCG" : "A",    "GCT" : "A",
    "GGA" : "G",    "GGC" : "G",    "GGG" : "G",    "GGT" : "G",
    "GTA" : "V",    "GTC" : "V",    "GTG" : "V",    "GTT" : "V",
    "TAA" : "-",    "TAC" : "Y",    "TAG" : "-",    "TAT" : "Y",
    "TCA" : "S",    "TCC" : "S",    "TCG" : "S",    "TCT" : "S",
    "TGA" : "-",    "TGC" : "C",    "TGG" : "W",    "TGT" : "C",
    "TTA" : "L",    "TTC" : "F",    "TTG" : "L",    "TTT" : "F",
    }

def one_frame_translate(seq, stub, id, this_frame):
    raw_orfs = []

    # Frame 0:
    frame = []
    for c in range(0, len(seq), 3):
        #cseq = seq[c:c+3]
        if 'N' in seq[c:c+3]:
            continue

        if len(seq[c:c+3]) != 3:
            break
        frame.append(codons[seq[c:c+3]])

    frame = ''.join(frame)

    #print(frame)

    orfs = []
    for m in re.finditer(r'\-[A-Z]*.\-', frame):
        start, item, end = m.start()+1, m.group(), m.end()-1
        if len(item) < min_length_to_consider:
            continue
        orfs.append((f"{stub};{id};{this_frame};{start*3}:{end*3}", item[1:-1]))

    #print(orfs)
    print(f"Processed id: {id} frame: {this_frame}, found {len(orfs):,}")
    return orfs

def three_frame_translation(seq, id, stub):
    raw_orfs = []
    raw_orfs += one_frame_translate(seq, stub, id, '0+')
    raw_orfs += one_frame_translate(seq[1:], stub, id, '1+')
    raw_orfs += one_frame_translate(seq[2:], stub, id, '2+')
    rev_comp = rc(seq)
    raw_orfs += one_frame_translate(rev_comp, stub, id, '0-')
    raw_orfs += one_frame_translate(rev_comp[1:], stub, id, '1-')
    raw_orfs += one_frame_translate(rev_comp[2:], stub, id, '2-')
    return raw_orfs

# Do all this the dumb way:

def process(infile, outfile, stub):
    print(f"Doing {infile}...")
    oh = gzip.open(infile, 'rt')

    orfs = []
    ids = []
    curr_id = None
    genome_seq = []
    pos = 0
    for line in oh:
        if line[0] == '>':
            # Process the last ID
            if curr_id:
                orfs += three_frame_translation(''.join(genome_seq), curr_id, stub)

            del genome_seq

            pos = 0 #Record current ID
            curr_id = line[1:].strip().upper().split(' ')[0] # hg38
            genome_seq = []
            print(f"Found ID '{curr_id}'")
            continue

        # do the line:
        line = line.strip()
        genome_seq.append(line)
    oh.close()

    # need to process the fall off:
    if curr_id:
        orfs += three_frame_translation(''.join(genome_seq), curr_id, stub)

    print(f"In total found: {len(orfs):,} ORFs")
    oh = gzip.open(outfile, 'wt')
    for orf in orfs:
        oh.write(f">{orf[0]}\n{orf[1]}\n")
    oh.close()
    return

if __name__ == '__main__':
    process('test.fa.gz', 'orfs_test.fa.gz', stub='test')
    process('Branchiostoma_lanceolatum.BraLan2.dna.toplevel.fa.gz', 'orfs_BraLan.fa.gz', stub='BraLan')
    process('Eptatretus_burgeri.Eburgeri_3.2.dna.toplevel.fa.gz', 'orfs_Eburgeri.fa.gz', stub='Eburgeri')
    process('Petromyzon_marinus.Pmarinus_7.0.dna.toplevel.fa.gz', 'orfs_Pmarinus.fa.gz', stub='Pmarinus')
    process('Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz', 'orfs_hg38.fa.gz', stub='hg38')

