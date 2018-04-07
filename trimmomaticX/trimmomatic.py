from statistics import mean
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class Trimmomatic:

    def __init__(self, path, out):
        self.path = ''
        self.out = ''
        self.sequences = []
        self.seqs = []

    def process_input(self):
        self.sequences = SeqIO.parse(self.path, 'fastq')

    def dump(self):
        SeqIO.write(self.seqs, self.out, 'fastq')

    def headcrop(self, sequence, size):
        # Drop size nucls from start
        # Create SeqRecord

        self.seqs.append(SeqRecord(Seq(str(sequence.seq)[size:]),
                                   id=sequence.id,
                                   name=sequence.name,
                                   description=sequence.description,
                                   letter_annotations={spec: vals[size:] for spec, vals in
                                                       sequence.letter_annotations.items()}))

        # sequence.letter_annotations['phred_quality'] = sequence.letter_annotations['phred_quality'][size:]
        # sequence.seq = str(sequence.seq)[size:]



    def tailcrop(self, sequence, size):
        # Drop size nucls from end

        self.seqs.append(SeqRecord(Seq(str(sequence.seq)[:-size]),
                                   id=sequence.id,
                                   name=sequence.name,
                                   description=sequence.description,
                                   letter_annotations={spec: vals[:-size] for spec, vals in
                                                       sequence.letter_annotations.items()}))

    def quality_crop(sequence, sliding_window, quality):
        #
        for i in range(len(sequence) - sliding_window + 1):
            if mean(sequence.letter_annotations['phred_quality'][i:i + sliding_window]) < quality:
                sequence.seq = str(sequence.seq)[:i]

    def pack_task(task, *args):
        return lambda x: task(x, *args)

    def do(sequences, tasks):
        for task in tasks:
            for sequence in sequences:
                task(sequence)