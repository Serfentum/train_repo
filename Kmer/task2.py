from fragmentator import Fragmentator
from Bio import SeqIO


with open('test.fasta', 'r') as file:
    data = SeqIO.parse(file, 'fasta')
    data = str(list(data)[0].seq)

a = Fragmentator(data)
a.sequencing(100).kmering(23)

a.kmer_plot(log_scale=True)
a.kmer_plot(60, 150)