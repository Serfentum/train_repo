from collections import Counter
from Bio import SeqIO
import matplotlib.pyplot as plt


class KmerSpectrum:

    def __init__(self, file_name, k=15, q=0):
        self.file_name = file_name
        self.k = k
        self.q = q
        self.kmers = Counter()
        self.data = None

    def analyze(self):
        sequences = SeqIO.parse(self.file_name, 'fastq')
        for record in sequences:
            self.fragmentate(record.seq, record.letter_annotations['phred_quality'])

    def fragmentate(self, read, quality):
        self.kmers.update(Counter([read[i:i + self.k] for i in range(len(read) - self.k + 1)
                                   if all(x >= self.q for x in quality[i:i + self.k])]))

    def transform(self):
        self.data = Counter(self.kmers.values())

    def genome_length_estimate(self):
        expectation = sum([i * j for i, j in zip(self.data.keys(), self.data.keys())])
        mp = self.maximum()
        print(f'Putative maximum is {mp}')
        return expectation / mp

    def maximum(self):
        """Try to find global maximum, but it could be local. If doesn`t find, take empiric value"""
        m = max(self.data.keys())
        for i in range(10, min(round(0.8 * m), m - 15)):
            for j in range(10):
                if self.data[i] < self.data[i - j] or self.data[i] < self.data[i + j]:
                    break
            else:
                return i
        return 0.4 * max(self.data.keys())


    def plot(self, xs=(None, None), ys=(None, None), log_scale=False):
        # Drawing
        plt.bar(self.data.keys(), self.data.values(),
                label=f'{"Logarithmic scale" if log_scale else "Linear scale"}\nk = {self.k}\nq = {self.q}')

        # Plot settings
        if log_scale:
            plt.yscale('log')
        plt.xlim(xs[0], xs[1])
        plt.ylim(ys[0], ys[1])

        # Plot decorations
        plt.title('Distribution of k-mers')
        plt.xlabel('Frequency of k-mers')
        plt.ylabel('# of distinct k-mers')
        plt.legend()

        # Showing
        plt.show()

path = '/home/arleg/data/test_kmer.fastq'
a = KmerSpectrum(path, k=11, q=20)
a.analyze()
a.transform()
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        a.plot()
print(a.genome_length_estimate())

b = KmerSpectrum(path, k=11, q=0)
b.analyze()
b.transform()
b.plot()
print(b.genome_length_estimate())



