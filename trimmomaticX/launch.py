#! /usr/bin/env python

from colored_argparse import ColoredArgParser
from trimmomatic import *



def parsing():
    # Parser initialization
    parser = ColoredArgParser(description='Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do '
                                          'eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim '
                                          'ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut '
                                          'aliquip ex ea commodo consequat.',
                              epilog='Goodbye',
                              fromfile_prefix_chars='@')
    parser.add_argument('-i', '--input', help='path to file', metavar='<path to file>',
                        type=str, required=True)
    parser.add_argument('-o', '--output', help='path to output', metavar='<path to output>',
                        type=str, default=name)
    parser.add_argument('-s', '--start', help='number of nucleotide to crop from start',
                        metavar='<1-sequence length>',
                        type=int)
    parser.add_argument('-e', '--end',
                        help='number of nucleotide to crop from end',
                        metavar='<1-sequence length>',
                        type=int)
    parser.add_argument('-w', '--sliding-window',
                        help='size of sliding window',
                        metavar='<1-sequence length>',
                        type=int)
    parser.add_argument('-q', '--quality',
                        help='necessary mean quality of fragment at sliding window',
                        type=int)
    parser.add_argument('-v', '--version', action='version', version='{} pre-alpha'.format(parser.prog))


    args = parser.parse_args()

    for arg, value in args.__dict__.items():
        print(type(arg), type(value), arg, value, sep='\t')

    print(args.__dict__)


    path = args.__dict__['input']
    out = args.__dict__['output']
    start = args.__dict__['start']
    end = args.__dict__['end']
    sliding_window = args.__dict__['sliding_window']
    quality = args.__dict__['quality']

    cleave = Trimmomatic(path, out)

    cleave.do([cleave.pack_task(cleave.headcrop, start),
               cleave.pack_task(cleave.tailcrop, end),
               cleave.pack_task(cleave.quality_crop, sliding_window, quality)])

    cleave.dump()


if __name__ == '__main__':
    name = 'trim_output.fastq'
    parsing()