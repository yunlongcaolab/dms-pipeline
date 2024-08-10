import sys
import pandas as pd
from random import random

use_barcodes = set(pd.read_csv(sys.argv[1])['barcode'].to_list())

# read fastq from stdin and filter

prob_missing = 0.01 # probability of including a mismatched barcode

line = sys.stdin.readline()
while line:
    if line[0] == '@':
        header = line
        seq = sys.stdin.readline()
        plus = sys.stdin.readline()
        qual = sys.stdin.readline()
        if seq.strip()[0:26] in use_barcodes or random() < prob_missing:
            sys.stdout.write(header)
            sys.stdout.write(seq)
            sys.stdout.write(plus)
            sys.stdout.write(qual)
    line = sys.stdin.readline()