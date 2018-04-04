from kmer_spectrum import *


for q in (0, 20):
    with open(f'{q}qual') as coords:
        xy = {}
        for line in coords:
            line = line.strip().split()
            xy[int(line[0])] = int(line[1])

        b = KmerSpectrum('', k=11, q=q)
        dat = xy
        b.data = {i: 0 for i in range(1, max(dat))}
        b.data.update(dat)

        for i in range(3):
            b.plot(f'plot_{q}_quality_{i}_smoothes.svg')
            print(f'Genome length estimate with {q} quality threshold and {i} smoothes - {b.genome_length_estimate()}')
            b.smooth_function()
