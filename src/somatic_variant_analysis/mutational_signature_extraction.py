from SigProfilerExtractor import sigpro as sig
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument(
    "-i", action="store", dest="input", help="Enter the filepath of your variants file"
)
parser.add_argument(
    "-o", action="store", dest="output", help="Enter the output file path"
)
parser.add_argument(
    "-g",
    action="store",
    dest="genome",
    help="Enter the genome used to align the sequencing data",
)
parser.add_argument(
    "-t",
    action="store",
    dest="threads",
    help="Enter the genome used to align the sequencing data",
    type=int,
)
args = parser.parse_args()

df = pd.read_csv(args.input, sep=",")

if __name__ == "__main__":
    sig.sigProfilerExtractor(
        "matrix",
        args.output,
        df,
        opportunity_genome=args.genome,
        exome=True,
        minimum_signatures=1,
        maximum_signatures=10,
        cpu=args.threads,
    )
