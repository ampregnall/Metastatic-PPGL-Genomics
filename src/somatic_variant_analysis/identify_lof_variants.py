"""
Created on Sun Nov 12 22:03:15 2023

@author: ampregnall
email: andrew.pregnall@pennmedicine.upenn.edu
"""

import csv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", action="store", dest="input", help="Enter the filepath of your variants file"
)
parser.add_argument(
    "-lof", action="store", dest="lof", help="Enter the reference genome"
)
parser.add_argument("-gof", action="store", dest="gof", help="Enter the project name")
parser.add_argument(
    "-af",
    action="store",
    dest="allele_freq",
    help="Enter the population allele frequency for filtering",
)
args = parser.parse_args()

with open(args.input) as infile:
    with open(args.lof, "w") as lof:
        with open(args.gof, "w") as gof:
            reader = csv.DictReader(infile, delimiter="\t")
            lof_writer = csv.DictWriter(
                lof, delimiter="\t", fieldnames=reader.fieldnames
            )
            gof_writer = csv.DictWriter(
                gof, delimiter="\t", fieldnames=reader.fieldnames
            )

            next(reader)
            lof_writer.writeheader()
            gof_writer.writeheader()

            for row in reader:
                if (
                    row["gnomAD.MAX_AF"] == "."
                    or float(row["gnomAD.MAX_AF"]) < float(args.allele_freq)
                ) and "benign" not in row["ClinVar.SIG"].lower():
                    if (
                        "frameshift_variant" in row["Variant.Consequence"]
                        or "stop_gained" in row["Variant.Consequence"]
                    ):
                        lof_writer.writerow(row)
                    elif "missense_variant" in row["Variant.Consequence"]:
                        if row["REVEL"] != ".":
                            if float(row["REVEL"]) > 0.5:
                                lof_writer.writerow(row)
                            else:
                                gof_writer.writerow(row)
                        else:
                            gof_writer.writerow(row)
                    elif (
                        row["Variant.Consequence"] == "splice_acceptor_variant"
                        or row["Variant.Consequence"] == "splice_donor_variant"
                    ):
                        lof_writer.writerow(row)
                    else:
                        gof_writer.writerow(row)
