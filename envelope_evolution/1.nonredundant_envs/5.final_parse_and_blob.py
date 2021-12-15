
# The aim here is to now output blobbed FASTAs for HMM generation.
# One remaining problem is that there is likely to still be trivial matches remaining in the search.

# TODO: Wouldn't it be better to shrink this to the smallest non-trivial list?

import sys, os, gzip
from glbase3.utils import convertFASTAtoDict

oh = gzip.open('simple_filtered.results.gz', "rt")

idx = 0

blobs = {}
per_search_blob = {}

for line in oh:
    idx += 1
    if idx % 1e5 == 0:
        print('Processed {:,}'.format(idx))
        #break

    l = line.strip().split(',')
    length = float(l[2]) # percent length matched

    #if length < 90: # 90% overlap
    #    # It's one of the trivial overlaps;
    #    continue

    if length < 60: # 60% overlap
        continue

    id1 = l[0]
    id2 = l[1]

    if id1 == id2: #Ignore 100% matches;
        continue

    if id1 in per_search_blob:
        blobs[per_search_blob[id1]].add(id2)
        per_search_blob[id2] = per_search_blob[id1]
    elif id2 in per_search_blob:
        blobs[per_search_blob[id2]].add(id1)
        per_search_blob[id1] = per_search_blob[id2]
    else:
        # Not in any blob, add a new one:
        blobs[idx] = set([id1, id2]) # key is not important and has no meaning
        per_search_blob[id1] = idx
        per_search_blob[id2] = idx
        # It's possible that a id can end up registered to two blobs,
        # But as I am going to do 2+ passes, seems not a big problem?

for k in blobs:
    print('Blob {} members: {}'.format(k, len(blobs[k])))
print('Number of blobs: {:,}'.format(len(blobs)))

# get all the fastas:
fasta = convertFASTAtoDict('simple_filtered.fasta')
fasta = {f['name'].split(' ')[0]: f['seq'] for f in fasta}
#print(fasta)

# save the ID of the first one in the blob;
oh = open('single_representative_envs.fasta', 'wt')
ids_saved = set([])
for k in blobs:
    # Find the k with the longest sequence
    longest_so_far = -1 # The fid
    longest_so_far_len = 0
    for fid in blobs[k]:
        if len(fasta[fid]) > longest_so_far_len:
            longest_so_far_len = len(fasta[fid])
            longest_so_far = fid

    if longest_so_far not in ids_saved:
        oh.write('>{}\n{}\n'.format(longest_so_far, fasta[longest_so_far]))
        ids_saved.add(longest_so_far)

oh.close()
