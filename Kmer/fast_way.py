from Bio import SeqIO


with open('seq_y_pestis.fasta', 'r') as file:
    res = {}
    data = SeqIO.parse(file, 'fasta')
    data = str(list(data)[0].seq)
    kmer_size = 23

    for i in range(len(data) - kmer_size + 1):
        res[data[i : i+kmer_size]] = res.get(data[i : i+kmer_size], 0) + 1

for i in sorted(res, key=lambda x: res[x], reverse=True):
    print(f'{i}\t{res[i]}')

print(f"\nYour answer is {sorted(res.items(), key=lambda x: x[1], reverse=True)[0]}")
