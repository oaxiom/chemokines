
# The output limit per search is 500 matches, so there is some significant banding
# and this will need to be done for at least 2 passes:

import sys, os, gzip
from glbase3.utils import convertFASTAtoDict

oh = gzip.open('all_vs_all.results.gz', "rt")

idx = 0

blobs = {}
per_search_blob = {}

for line in oh:
    idx += 1
    if idx % 1e6 == 0:
        print('Processed {:,}'.format(idx))
        #break

    l = line.strip().split(',')
    length = float(l[2])

    if length < 90:
        continue

    id1 = l[0]
    id2 = l[1]

    '''
    # Comprehensive, but too slow past >48,000,000 matches;
    # see if it's already in a blob, and then add it to that set
    for blob in blobs:
        if id1 in blobs[blob]:
            blobs[blob].add(id2)
            break
        if id2 in blobs[blob]:
            blobs[blob].add(id1)
            break
    else:
        # Not in any blob, add a new one:
        blobs[idx] = set([id1, id2]) # key is not important and has no meaning
    '''

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
fasta = convertFASTAtoDict('1.blastdb/all_envs.fasta')
fasta = {f['name'].split(' ')[0]: f['seq'] for f in fasta}
#print(fasta)

# save the ID of the first one in the blob;
oh = open('simple_filtered.fasta', 'wt')
for k in blobs:
    # Find the k with the longest sequence
    longest_so_far = -1 # The fid
    longest_so_far_len = 0
    for fid in blobs[k]:
        if len(fasta[fid]) > longest_so_far_len:
            longest_so_far_len = len(fasta[fid])
            longest_so_far = fid

    oh.write('>{}\n{}\n'.format(longest_so_far, fasta[longest_so_far]))
oh.close()
