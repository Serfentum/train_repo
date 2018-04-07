#! /usr/bin/env python

from colored_argparse import ColoredArgParser
from trimmomaticx import *


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
    parser.add_argument('-e', '--end', help='number of nucleotide to crop from end',
                        metavar='<1-sequence length>')
    parser.add_argument('-w', '--sliding-window', help='another flag', action='store_false')
    parser.add_argument('-q', '--quality', help='another flag', action='store_false')
    parser.add_argument('-v', '--version', action='version', version='{} pre-alpha'.format(parser.prog))


    args = parser.parse_args()

    for arg, value in args.__dict__.items():
        print(type(arg), type(value), arg, value, sep='\t')

        # a = process_input(args.__dict__['input'])
        #
        # t = pack_task(quality_crop, 5, 3)



if __name__ == '__main__':
    name = 'trim_output.fastq'
    parsing()