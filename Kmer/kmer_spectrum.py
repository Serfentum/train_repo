import numpy as np
from collections import Counter
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
        self.derivative = None
        self.threshold = None
        self.max = None

    def analyze(self):
        """Read fastq file, fragmentate it into kmers, take only kmers with all bases with certain quality"""
        sequences = SeqIO.parse(self.file_name, 'fastq')
        for record in sequences:
            self.fragmentate(record.seq, record.letter_annotations['phred_quality'])

    def fragmentate(self, read, quality):
        """
        Fragmentate read to kmers, add kmers to Counter() if the quality of all letters >= quality threshold
        :param read: str - read sequence
        :param quality: list - with ints for every base
        :return:
        """
        self.kmers.update(Counter([read[i:i + self.k] for i in range(len(read) - self.k + 1)
                                   if all(x >= self.q for x in quality[i:i + self.k])]))

    def transform(self):
        """Count number of different reads with same occurency, make y of missing x in dict equal to 0"""
        information = Counter(self.kmers.values())

        # Create np array of plot coordinates from dict
        distribution = np.zeros(max(information) + 1)
        for position, value in information.items():
            distribution[position] = value
        self.data = distribution

    def genome_length_estimate(self):
        """
        Estimate length of genome upon kmer distribution and estimation of main peak abscissa
        :return: float - length of genome in bases (perhaps I should round it to int)
        """
        if not self.threshold:
            self.cutoff()
        self.maximum()
        # Compute sum of x * y for all x and divide by x of main peak and 2 due 2 double strandness of coming DNA
        expectation = np.sum(np.arange(self.threshold, len(self.data)) * self.data[self.threshold:]) / 2
        return expectation / self.max

    def unit_derivative(self):
        """
        Here I compute derivative of function with increment of x equal to 1
        :return:
        """
        # Compute f(x) - f(x - 1) - derivative of y
        self.derivative = self.data[2:] - self.data[1:-1]

    def cutoff(self):
        """
        # Find x by which we will separate noise
        :return:
        """
        # Compute derivative
        self.unit_derivative()

        # Find minimum - first point where function start grow after stable period
        # #TODO perhaps correct function, though looks ok according to tests
        self.threshold = np.where(self.derivative > 0)[0][0] - 1

    def maximum(self):
        """
        Find x coordinate of main peak
        :return:
        """
        # Find position of main peak
        self.max = np.where(self.derivative[self.threshold:] < 0)[0][0] + self.threshold


    def smooth_function(self, window=10):
        """
        Smooth our function of number of distinct kmers depending upon frequency of kmers.
        Slides on x-axis and substitute y-values on median at window.
        Update self.data
        :param window: int - length of sliding window
        :return:
        """
        smoothed = np.zeros(len(self.data))
        for i in range(len(self.data)):
            smoothed[i] = np.median(self.data[i:i + window])
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
        plt.bar(range(1, len(self.data) + 1), self.data,
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
    path = '/home/arleg/data/ttest_kmer.fastq'

    for q in (0,):
        b = KmerSpectrum(path, k=11, q=q)
        b.analyze()
        b.transform()
        for i in range(2):
            b.smooth_function()
            b.plot(f'NewPlot_{q}_quality_{i}_smoothes.png')
            print(f'Genome length estimate with {q} quality threshold and {i} smoothes - {b.genome_length_estimate()}')
            print(b.threshold)
            print(b.max)



