#! /usr/bin/env python

from colored_argparse import ColoredArgParser
from trimmomatic import *



def run():
    # Parser initialization
    parser = ColoredArgParser(description='Tool for cleaving fastq reads',
                              epilog='Goodbye')
    parser.add_argument('-i', '--input',
                        help='path to file',
                        metavar='<path to file>',
                        type=str,
                        required=True)
    parser.add_argument('-o', '--output',
                        help='path to output',
                        metavar='<path to output>',
                        type=str,
                        default=name)
    parser.add_argument('-s', '--start',
                        help='number of nucleotide to crop from start',
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
                        help='threshold of mean fragment quality at sliding window',
                        metavar='<phred quality>',
                        type=int)
    parser.add_argument('-v', '--version',
                        action='version',
                        version='{} pre-alpha'.format(parser.prog))

    args = parser.parse_args()
    # Make options aliases
    path = args.__dict__['input']
    out = args.__dict__['output']
    start = args.__dict__['start']
    end = args.__dict__['end']
    sliding_window = args.__dict__['sliding_window']
    quality = args.__dict__['quality']

    # Initialize class
    cleave = Trimmomatic(path, out)
    # Apply all possible functions with appropriate arguments to sequences
    cleave.do([cleave.pack_task(cleave.headcrop, start),
               cleave.pack_task(cleave.tailcrop, end),
               cleave.pack_task(cleave.quality_crop, sliding_window, quality)])
    # Write to file
    cleave.dump()

if __name__ == '__main__':
    # Give default name for output and run processing
    name = 'trim_output.fastq'
    run()