#! /usr/bin/env python

from colored_argparse import ColoredArgParser


def parsing():
    # Parser initialization
    parser = ColoredArgParser(description='Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do '
                                          'eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim '
                                          'ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut '
                                          'aliquip ex ea commodo consequat.',
                              epilog='Goodbye',
                              fromfile_prefix_chars='@')
    parser.add_argument('-i', '--input', help='absolute or relative path to file', metavar='<path to file>',
                        type=str, required=True)
    parser.add_argument('-o', '--output', help='output description', metavar='<path to output>',
                        type=str, default=name)
    parser.add_argument('-b', '--binary', help='some flag', action='store_true')
    parser.add_argument('-t', '--thetary', help='numerical argument', metavar='<1-10>')
    parser.add_argument('-z', '--zetary', help='another flag', action='store_false')
    parser.add_argument('-v', '--version', action='version', version='{} pre-alpha'.format(parser.prog))


    args = parser.parse_args()

    for arg, value in args.__dict__.items():
        print(arg, value, sep='\t')


if __name__ == '__main__':
    name = 'trim_output.fastq'
    parsing()