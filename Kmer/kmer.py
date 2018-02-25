class Kmer:
    """
    Class for kmer. Instance for each kmer should be instantiated.
    I didn`t make setters and getters cause it`s more common practice to get attribute value directly in python.
    """
    sequence = ''
    occurence = 0
    start_positions = []

    def __init__(self, sequence='', occurence=0, start_positions=[]):
        self.sequence = sequence
        self.occurence = occurence
        self.start_positions = start_positions

    def _occured(self, n=1):
        self.occurence += n

    def _add_locus(self, position):
        self.start_positions.append(position)

    def _find_kmer(self, seq, start=None, stop=None):
        pos = seq.find(self.sequence, start, stop)
        if pos != -1:
            self._add_locus(pos)
            self._occured()
        return pos

    def find_all_kmer_occurences(self, seq, start=None, stop=None):
        """Find all occurences of kmer in sequence from start index to stop. Update list with start positions.
        Sequential usage of this method on 1 object with intersected search regions could led to redundant list.
        Yet I can write getter to start_positions and remove duplicates on call"""
        perhaps_part_of_positions = []
        pos = self._find_kmer(seq, start, stop)
        while pos != -1:
            perhaps_part_of_positions.append(pos)
            pos = self._find_kmer(seq, pos + 1, stop)
        return perhaps_part_of_positions



    def dump(self):
        """
        Print tsv with kmer sequence and all its position where it placed.
        """
        for pos in sorted(self.start_positions):
            print(f'{self.sequence}\t{pos}')


# Test. We should write unit tests
# s = 'ATGGGGCCGGTG'
# a = Kmer('GGG')
# print(a.find_all_kmer_occurences(s))
# a.dump()
# print(a.occurence)
# print(a.start_positions)