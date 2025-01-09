"""
Created on Sun Nov 12 18:46:15 2023

@author: ampregnall
email: andrew.pregnall@pennmedicine.upenn.edu
"""

import csv
import os
import argparse

# Take in a list of variants, which can be one tumor's unfiltered file or many tumors' variants. Must be concatenated into one file.
parser = argparse.ArgumentParser()

parser.add_argument(
    "-i", action="store", dest="input", help="Enter the filepath of your variants file"
)
parser.add_argument(
    "-o", action="store", dest="output", help="Enter the output file path"
)
parser.add_argument(
    "-rd", action="store", dest="read_depth", help="Enter the read depth for filtering"
)

args = parser.parse_args()

# Check input
if args.input == None:
    print("Error: please specify an input file of variants with -i")
    exit()
if not os.path.exists(args.input):
    print("Input file does not exist")
    exit()

# Read each line of input into the output based on the desired filtering criteria. Exclude rows that do not meet all filtering criteria.
with open(args.input, "r") as infile, open(args.output, "w") as outfile:
    reader = csv.DictReader(infile, delimiter="\t")
    writer = csv.DictWriter(
        outfile, delimiter="\t", fieldnames=reader.fieldnames, lineterminator="\r"
    )
    writer.writeheader()

    for row in reader:
        if float(row["Tumor.AltDepth"]) >= float(args.read_depth):
            writer.writerow(row)
