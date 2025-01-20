import os
import argparse
import csv

# Take in a list of variants, which can be one tumor's unfiltered file or many tumors' variants. Must be concatenated into one file.
parser = argparse.ArgumentParser()
parser.add_argument('-i', action='store', dest='input', help='Enter the filepath of your variants file')
parser.add_argument("-o", action='store', dest="output", help='Enter the output file path')
parser.add_argument("-p", action="store", dest='project', help='Enter the project name')
parser.add_argument("-rd", action='store', dest='read_depth', help='Enter the read depth for filtering')
parser.add_argument("-g", action="store", dest='genome', help='Enter the genome used to align the sequencing data')
args = parser.parse_args()

# Check input
if args.input == None:
    print("Error: please specify an input file of variants with -i")
    exit()
if not os.path.exists(args.input):
    print("Input file does not exist")
    exit()

# Construct a header for the output file; only sample (1), chrom (5), start (6), end (7), ref (8), and mut (9) are required
outfields = ['Project', 'Sample', 'ID', 'Genome', 'mut_type', 'chrom', 'pos_start', 'pos_end', 'ref', 'alt', 'type']

# Read each line of input into the output based on the desired filtering criteria. Exclude rows that do not meet all filtering criteria.
with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
    reader = csv.DictReader(infile, delimiter=',')
    writer = csv.DictWriter(outfile, delimiter='\t', fieldnames=outfields, lineterminator='\r')
    writer.writeheader()

    for row in reader:
        
        if float(row['Tumor.AltDepth']) >= float(args.read_depth):

            newDict = {'Project': args.project,
                    'Sample': row['Tumor.ID'],
                    'ID': '.', # Assume this is a unique variant ID? Not needed for now
                    'Genome': args.genome,
                    'mut_type': row['Variant.Class'],
                    'chrom': row['Chr'],
                    'pos_start': row['Start'],
                    'pos_end': row['Start'], # Don't think this matters based on input
                    'ref': row['REF'],
                    'alt': row['ALT'],
                    'type': 'SOMATIC'}
            
            writer.writerow(newDict)
