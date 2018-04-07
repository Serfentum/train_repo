from statistics import mean
from Bio import SeqIO


class Trimmomatic:

    def __init__(self, path, out='~/train_repo/trimmomaticX/trimmed.fastq'):
        self.sequences = list(SeqIO.parse(path, 'fastq'))
        self.out = out
        self.seqs = []

    def dump(self):
        SeqIO.write(self.sequences, self.out, 'fastq')

    def headcrop(self, sequence, size):
        self.seqs.append(sequence[size:])

    def tailcrop(self, sequence, size):
        self.seqs.append(sequence[:-size])

    def quality_crop(self, sequence, sliding_window, quality):
        for i in range(len(sequence) - sliding_window + 1):
            if mean(sequence.letter_annotations['phred_quality'][i:i + sliding_window]) < quality:
                break
        self.seqs.append(sequence[0:i + sliding_window]
                         if mean(sequence.letter_annotations['phred_quality'][i:i + sliding_window]) >= quality
                         else sequence[0:i])

    def pack_task(self, task, *args):
        if [i for i in args if i]:
            return lambda x: task(x, *args)

    def do(self, tasks):
        for task in tasks:
            if task:
                for sequence in self.sequences:
                    task(sequence)

                self.sequences = self.seqs
                self.seqs = []


if __name__ == '__main__':
    a = Trimmomatic('/home/arleg/train_repo/trimmomaticX/test', 'out')

    a.do([a.pack_task(a.tailcrop, 5), a.pack_task(a.headcrop, 5), a.pack_task(a.quality_crop, 3, 38)])

    a.dump()
