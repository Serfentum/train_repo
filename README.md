# Genebank parser and Kmer class
This tool is devoted to parse output from genebank and was made as a homework to your assignment, Jenya. You can find all records which contains one of search words in your query.

## Prerequisites
To run our tool you will need
- Python Interpretator 3.X version
- argparse library

## Usage
`$ python parse_gb.py options`

Options have 2 forms - short 1-letter and long. Long is prefixed with `--` and short with `-` Possible options (in form short long - description)
- h help	- display help information about tool
- i input	- phrase after it should be path to your input file from genbank, which you`ll parse
- o output	- phrase after it should be path to output file with its name
- s specs	- words which will be searched in input
- v version	- display information about version of tool
After writing options you should press button `Enter` and wait until searching will be done. Output file will be written in specified location. By default it is current directory and file with name *parsedGB.csv*

## Authors
Arleg

## License
BI License

## Acknowledgements
I`d like to thanks everybody
