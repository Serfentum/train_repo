from collections import Counter
import matplotlib.pyplot as plt


class Fragmentator:
    """Class for sequencing and fragmenting kmers"""

    def __init__(self, seq):
        self.seq = seq
        self.reads = []
        self.kmers = []
        self.points = {}

    def sequencing(self, read_length=100):
        """Sequence sequence. Coverage of ends is less than other part has"""
        self.reads = [self.seq[i:i + read_length] for i in range(len(self.seq) - read_length + 1)]
        return self

    def kmering(self, k=23):
        """Fragment kmers from reads"""
        self.k = k
        self.kmers = [read[i:i + k] for read in self.reads for i in range(len(read) - k + 1)]
        return self

    def kmer_stats(self):
        """Count coordinates of points for kmer plot"""
        # Frequency of kmers
        kmer_dict = dict(Counter(self.kmers))

        # Frequency of kmer frequencies
        self.points = dict(Counter(kmer_dict.values()))
        return self


    def kmer_plot(self, mi=None, ma=None, log_scale=False):
        """Plot kmer frequencies distribution"""
        # Compute coordinates if we don`t have them
        if not self.points:
            self.kmer_stats()

        # Plotting
        plt.bar(self.points.keys(), self.points.values(), label=f'k = {self.k}')
        # Zooming
        plt.xlim(mi, ma)
        plt.title('Distribution of k-mers')
        plt.xlabel('Frequency of k-mers')
        plt.ylabel('# of distinct k-mers')
        if log_scale:
            plt.yscale('log')
        plt.legend()
        plt.show()


# a = Fragmentator('ATGCCTACCTGGTCCTACGCCT')
# a.sequencing(5).kmering(3).kmer_plot()
# a.sequencing(5).kmering(3).kmer_plot(2)
# a.sequencing(5).kmering(3).kmer_plot(ma=5)
# a.sequencing(5).kmering(3).kmer_plot(mi=3, ma=5)