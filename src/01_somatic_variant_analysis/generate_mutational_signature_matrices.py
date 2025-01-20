#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 19 16:10:09 2023

@author: ampregnall
email: andrew.pregnall@pennmedicine.upenn.edu
"""

import argparse
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

parser = argparse.ArgumentParser()
parser.add_argument('-i', action='store', dest='input', help='Enter the filepath of your variants file')
parser.add_argument("-g", action='store', dest="genome", help='Enter the reference genome')
parser.add_argument("-p", action='store', dest='project', help="Enter the project name")
args = parser.parse_args()

matrices = matGen.SigProfilerMatrixGeneratorFunc(args.project, args.genome, args.input,
                                                 plot=True, exome=False, bed_file=None, chrom_based=False, 
                                                 tsb_stat=True, seqInfo=True, cushion=100)