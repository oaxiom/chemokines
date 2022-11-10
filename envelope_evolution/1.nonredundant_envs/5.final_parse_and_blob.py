
# The aim here is to now output blobbed FASTAs for HMM generation.
# One remaining problem is that there is likely to still be trivial matches remaining in the search.

# TODO: Wouldn't it be better to shrink this to the smallest non-trivial list?

# The assumption here is that every env has a hit.
# If it has no hits then it wont appear in the simple_filtered file?

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

    if length < 60: # at least x0% of the sequence overlaps;
        # Higher will bundle more together, lower will split the envs up and produce more candidates;
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

for k in blobs:
    print('Blob {} members: {}'.format(k, len(blobs[k])))
print('Number of blobs: {:,}'.format(len(blobs)))

# get all the fastas:
fasta = convertFASTAtoDict('simple_filtered.fasta')
fasta = {f['name'].split(' ')[0]: f['seq'] for f in fasta}
#print(fasta)

# Make certain all envs appear at least once anywhere in the results.
# If it's not found, add it as an extra fasta entry to save.

_found = 0
_not_found = 0
toadd = []
for f in fasta:
    # Is it in anyblob?
    for k in blobs:
        if f in blobs[k]:
            _found += 1
            break
    else:
        _not_found += 1
        toadd.append(f)

print(f'Number of envs not found: {_not_found} presumed unique')

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

# And add all ths missing (unique?) envs:
for f in toadd:
    oh.write('>{}\n{}\n'.format(f, fasta[f]))

oh.close()
