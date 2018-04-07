from collections import Counter
from statistics import median
from Bio import SeqIO
import matplotlib.pyplot as plt


class KmerSpectrum:
    """Class for analyzing kmer spectrum of sequence experiment"""
    def __init__(self, file_name, k=15, q=0):
        """
        Constructor
        :param file_name: str - path to fastq file
        :param k: int - kmer length
        :param q: int - quality threshold
        """
        self.file_name = file_name
        self.k = k
        self.q = q
        self.kmers = Counter()
        self.data = None
        self.derivative = {}
        self.threshold = None
        self.max = None

    def analyze(self):
        """Read fastq file, fragmentate it into kmers, take only kmers with all bases with certain quality"""
        sequences = SeqIO.parse(self.file_name, 'fastq')
        for record in sequences:
            self.fragmentate(record.seq, record.letter_annotations['phred_quality'])

    def fragmentate(self, read, quality):
        """
        Fragmentate read to kmers, add kmers to Counter() if the quality of all kmers >= quality threshold
        :param read: str - read sequence
        :param quality: list - with ints for every base
        :return:
        """
        self.kmers.update(Counter([read[i:i + self.k] for i in range(len(read) - self.k + 1)
                                   if all(x >= self.q for x in quality[i:i + self.k])]))

    def transform(self):
        """Count number of different reads with same occurency, make y of missing x in dict equal to 0"""
        information = Counter(self.kmers.values())
        self.data = {i: 0 for i in range(max(information))}
        self.data.update(information)

    def genome_length_estimate(self):
        """
        Estimate length of genome upon kmer distribution and estimation of main peak abscissa
        :return: float - length of genome in bases (perhaps I should round it to int)
        """
        self.maximum()
        expectation = sum([i * j for i, j in zip(self.data.keys(), self.data.values())])
        return expectation / self.max

    def unit_derivative(self):
        """
        Here I compute derivative of function with increment of x equal to 1
        :return:
        """
        for i in range(2, max(self.data) + 1):
            self.derivative[i] = self.data[i] - self.data[i - 1]

    def cutoff(self):
        """
        # Find x by which we will separate noise
        :return:
        """
        # Compute derivative
        self.unit_derivative()
        # Find minimum
        for i in range(2, len(self.derivative)):
            if self.derivative[i] == 0 and self.derivative[i + 1] > 0:
                self.threshold = i
                break

    def maximum(self):
        """
        Find x coordinate of main peak
        :return:
        """
        # Find position of main peak
        # self.max = max(self.derivative, key=lambda x: self.derivative[x])
        for i in range(self.threshold, len(self.derivative)):
            if self.derivative[i] >= 0 and self.derivative[i + 1] < 0:
                print(f'{self.q} - Putative maximum is {self.max}')
                self.max = i
                break

    def smooth_function(self, window=10):
        """
        Smooth our function of number of distinct kmers depending upon frequency of kmers.
        Slides on x-axis and substitute y-values on median at window.
        Update self.data
        :param window: int - length of sliding window
        :return:
        """
        smoothed = {}
        for i in range(1, max(self.data) - window + 1):
            smoothed[i] = median(self.data[j] for j in range(i, i + window))
        self.data = smoothed

    def plot(self, output, xs=(None, None), ys=(None, None), log_scale=False):
        """
        Plot kmer spectrum
        :param output: str - path to output file
        :param xs: tuple - range of x-axis
        :param ys: tuple - range of y-axis
        :param log_scale: boolean - whether y-axis is logarithmically scaled
        :return:
        """
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

        # Plot noise separator
        self.cutoff()
        plt.axvline(self.threshold, color='red')

        # Save plot
        plt.savefig(output)
        plt.close()


if __name__ == '__main__':
    path = 'test.fastq'

    for q in (0, 20):
        b = KmerSpectrum(path, k=11, q=q)
        b.analyze()
        b.transform()
        for i in range(2):
            b.smooth_function()
        b.plot('ttt')

        # for i in range(3):
        #     b.plot(f'TEST_plot_{q}_quality_{i}_smoothes.png')
        #     print(f'Genome length estimate with {q} quality threshold and {i} smoothes - {b.genome_length_estimate()}')
        #     b.smooth_function()




