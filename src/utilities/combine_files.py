"""
Created on Thu May 11 12:35:15 2023

@author: ampregnall
email: andrew.pregnall@pennmedicine.upenn.edu
"""

import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", action="store", dest="input", help="Enter the filepath of your variants file"
)
parser.add_argument(
    "-o", action="store", dest="output", help="Enter the output file path"
)
args = parser.parse_args()

filenames = os.listdir(args.input)
combined = pd.concat(
    [
        pd.read_csv(
            "{0}/{1}".format(args.input, f),
            sep=",",
            encoding="ISO-8859-1",
            low_memory=False,
        )
        for f in filenames
    ],
    ignore_index=True,
)

combined.drop(combined.filter(regex="Unname"), axis=1, inplace=True)
combined.to_csv(args.output, index=False, sep="\t")
