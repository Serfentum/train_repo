from Bio import SeqIO
from kmer import *

# Dictionaries to store kmers and occurences in
res = {}
rr = {}

# Read sequence from fasta
with open('seq_y_pestis.fasta', 'r') as file:
    data = SeqIO.parse(file, 'fasta')
    data = str(list(data)[0].seq)
    kmer_size = 23

    # Create kmer from each new fragment, count it occurence
    for i in range(len(data) - kmer_size + 1):
        if data[i:i + kmer_size] not in res:
            km = Kmer(data[i:i + kmer_size])
            res[data[i:i + kmer_size]] = km
            km.find_all_kmer_occurences(data)

    # Add kmer occurences in 2nd dictionary
    for seq, km in res.items():
        rr[seq] = km.occurence

# Print answer
for r in sorted(rr, key=lambda x: rr[x], reverse=True):
    print(f'{r}\t{rr[r]}')

print(f"\nYour answer is {sorted(rr.items(), key=lambda x: x[1], reverse=True)[0]}")
