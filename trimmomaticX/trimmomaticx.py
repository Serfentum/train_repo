from statistics import mean
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# input - fastq
# headcrop - slice from start
# tailcrop - slice from end
# sliding window - size of window
# quality - slice if mean quality of several nucleotides drops

seqs = []

def process_input(path):
    sequences = SeqIO.parse(path, 'fastq')
    return sequences

def dump(sequences, path):
    SeqIO.write(sequences, path, 'fastq')

def headcrop(sequence, size):
    # Drop size nucls from start
    # Create SeqRecord
    seqs.append(sequence[size:])
    # a = SeqRecord(Seq(str(sequence.seq)[size:]),
    #               id=sequence.id,
    #               name=sequence.name,
    #               description=sequence.description,
    #               letter_annotations={spec: vals[size:] for spec, vals in sequence.letter_annotations.items()})

    # sequence.letter_annotations['phred_quality'] = sequence.letter_annotations['phred_quality'][size:]
    # sequence.seq = str(sequence.seq)[size:]
    # seqs.append(a)



def tailcrop(sequence, size):
    # Drop size nucls from end
    seqs.append(sequence[:-size])
    # a = SeqRecord(Seq(str(sequence.seq)[:-size]),
    #               id=sequence.id,
    #               name=sequence.name,
    #               description=sequence.description,
    #               letter_annotations={spec: vals[:-size] for spec, vals in sequence.letter_annotations.items()})
    # seqs.append(a)

def quality_crop(sequence, sliding_window, quality):
    #
    for i in range(len(sequence) - sliding_window + 1):
        if mean(sequence.letter_annotations['phred_quality'][i:i + sliding_window]) < quality:
            j = i
            break
    seqs.append(sequence[:i])

    # a = SeqRecord(Seq(str(sequence.seq)[:i]),
    #               id=sequence.id,
    #               name=sequence.name,
    #               description=sequence.description,
    #               letter_annotations={spec: vals[:i] for spec, vals in sequence.letter_annotations.items()})
    # seqs.append(a)



def pack_task(task, *args):
    return lambda x: task(x, *args)

def do(sequences, tasks):
    for task in tasks:
        for sequence in sequences:
            task(sequence)


# a = process_input('test')
# t = pack_task(quality_crop, 5, 30)
# do(a, [t])
# for i in a:
#     print(i)
# for i, j in zip(SeqIO.parse('test', 'fastq'), seqs):
#     print(i, j, sep='\t')
#     # for x, y in zip(i, j):
#     #     print(x == y, x, y, sep='\t')
# print(seqs)
