from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt


def gc_content(path='test_classwork2.fastq',
               out_path='gc_cont',
               threshold=20):
    """
    Plot distribution of GC content in reads from nucleotides with appropriate quality or from  all nucleotides
    in fasta (equivalent to threshold == 0). Also prints average GC%.
    :param path: str - path to input fasta or fastq file
    :param out_path: str - name of output plot
    :param threshold: int - value of threshold which should be reached by nucleotide quality to count it
    :return:
    """
    gc_contents = []
    # Infer data format, load it
    # For each read compute GC% from all nucleotides or from nucleotides with quality higher than threshold
    with open(path, 'r') as source:
        if source.readline().startswith('>'):
            reads = SeqIO.parse(path, 'fasta')
            for read in reads:
                gc_contents.append(read.seq.count('G') + read.seq.count('C') / len(read.seq))

        else:
            reads = SeqIO.parse(path, 'fastq')
            for read in reads:
                qual = np.array(read.letter_annotations['phred_quality'])
                inds = np.where(qual >= threshold)[0]

                letters = [read.seq[i] for i in inds]
                gc_contents.append(letters.count('G') + letters.count('C') / len(letters))

    print(np.average(gc_contents))

    # Plot distribution of GC% in reads
    plt.hist(gc_contents)
    plt.title(f'GC content of reads, {threshold}')
    plt.xlabel('GC%')
    plt.ylabel('Occurency')
    plt.savefig(out_path)


if __name__ == '__main__':
    gc_content()
