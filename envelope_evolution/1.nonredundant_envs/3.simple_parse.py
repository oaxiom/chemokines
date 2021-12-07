
# The output limit per search is 500 matches, so there is some significant banding
# and this will need to be done for at least 2 passes:

import sys, os, gzip
from glbase3.utils import convertFASTAtoDict

oh = gzip.open('all_vs_all.results.gz', "rt")

idx = 0

blobs = {}

for line in oh:
    idx += 1
    if idx % 1e7 == 0:
        print('Processed {:,}'.format(idx))
        break

    l = line.strip().split(',')
    length = float(l[2])

    if length < 90:
        continue

    id1 = l[0]
    id2 = l[1]

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


for k in blobs:
    print('Blob {} members: {}'.format(k, len(blobs[k])))
print('Number of blobs: {:,}'.format(len(blobs)))

# get all the fastas:
fasta = convertFASTAtoDict('blastdb/all_envs.fasta')
fasta = {f['name'].split(' ')[0]: f['seq'] for f in fasta}
#print(fasta)

# save the ID of the first one in the blob;
oh = open('simple_filtered.fasta', 'wt')
for k in blobs:
    fid = blobs[k].pop()
    oh.write('>{}\n{}\n'.format(fid, fasta[fid]))
oh.close()
