from fragmentator import Fragmentator
from Bio import SeqIO


with open('seq_y_pestis.fasta', 'r') as file:
    data = SeqIO.parse(file, 'fasta')
    data = str(list(data)[0].seq)

a = Fragmentator(data)
# Sequencing with read length = 100 is very long
a.sequencing(len(data)).kmering(23)

a.kmer_plot(log_scale=True)
