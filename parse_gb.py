import re
import csv
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@',
                                     epilog='Bye-bye!)')
    parser.add_argument('-i', '--input', help='absolute or relative path to file with its name',
                        metavar='<path to file>', required=True)
    parser.add_argument('-o', '--output', help='absolute or relative path to output with its name',
                        metavar='<path to file>', default='parsedGB.csv')
    parser.add_argument('-s', '--specs', help='tags information about which you wanna find in GenBank file',
                        metavar='<options>', nargs='+')
    parser.add_argument('-v', '--version', action='version', version='{} version is 0.03'.format(parser.prog))

    # Parse it out
    args = parser.parse_args()
    path_to = args.input
    path_out = args.output
    tags = args.specs

    # Create choices
    tags = '|'.join(tags)

    # Fields to extract, linked with regexp
    fields = ('locus', 'def', 'accession', 'source', 'title')

    # Pattern for information finding
    pat = re.compile(r'''LOCUS(?P<locus>.*?)         # locus
                        DEFINITION(?P<def>.*?)       # definition
                        ACCESSION(?P<accession>.*?)  # accession
                        VERSION                      # skip version
                        .*?
                        SOURCE(?P<source>.*?)        # source
                        ORGANISM                     # skip organism
                        .*?
                        TITLE(?P<title>[^(JOURNAL)]*?# title
                        ''' + '({})'.format(tags) +  # insert tags to find
                        '.*?)JOURNAL',               # up to journal
                        re.X | re.DOTALL)

    # Pattern for whitespaces substitution
    strip = re.compile(r'''\s+''')

    # Read`n`Write
    with open(path_to) as file, open(path_out, 'w') as dest:
        file = file.read()
        it = re.finditer(pat, file)
        data = []

        # Collect fields, clean them with auxiliary pattern
        for x in it:
            data.append([re.sub(strip, ' ', x.group(i).strip()) for i in fields])

        # Write data to csv file
        dest = csv.writer(dest, quoting = csv.QUOTE_ALL)
        dest.writerow(fields)
        for i in data:
            dest.writerow(i)




