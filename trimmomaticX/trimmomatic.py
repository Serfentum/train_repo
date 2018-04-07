from statistics import mean
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class Trimmomatic:

    def __init__(self, path, out='~/train_repo/trimmomaticX/trimmed.fastq'):
        self.sequences = list(SeqIO.parse(path, 'fastq'))
        self.out = out
        self.seqs = []

    def dump(self):
        SeqIO.write(self.sequences, self.out, 'fastq')

    def headcrop(self, sequence, size):
        # Drop size nucls from start
        # Create SeqRecord
        self.seqs.append(sequence[size:])
        # self.seqs.append(SeqRecord(Seq(str(sequence.seq)[size:]),
        #                            id=sequence.id,
        #                            name=sequence.name,
        #                            description=sequence.description,
        #                            letter_annotations={spec: vals[size:] for spec, vals in
        #                                                sequence.letter_annotations.items()}))

        # sequence.letter_annotations['phred_quality'] = sequence.letter_annotations['phred_quality'][size:]
        # sequence.seq = str(sequence.seq)[size:]

    def tailcrop(self, sequence, size):
        # Drop size nucls from end
        self.seqs.append(sequence[:-size])
        # self.seqs.append(SeqRecord(Seq(str(sequence.seq)[:-size]),
        #                            id=sequence.id,
        #                            name=sequence.name,
        #                            description=sequence.description,
        #                            letter_annotations={spec: vals[:-size] for spec, vals in
        #                                                sequence.letter_annotations.items()}))

    def quality_crop(self, sequence, sliding_window, quality):
        for i in range(len(sequence) - sliding_window + 1):
            print(sequence.letter_annotations['phred_quality'][i:i + sliding_window])
            if mean(sequence.letter_annotations['phred_quality'][i:i + sliding_window]) < quality:
                j = i
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


a = Trimmomatic('/home/arleg/train_repo/trimmomaticX/test', 'out')
g = a.sequences # seq

for i in a.sequences:
    print(i.letter_annotations['phred_quality'])

print(list(g))

a.do([a.pack_task(a.tailcrop, 5), a.pack_task(a.headcrop, 5), a.pack_task(a.quality_crop, 3, 38)])
print(a.sequences)
a.dump()
