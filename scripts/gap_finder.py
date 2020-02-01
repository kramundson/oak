import sys
import os
from re import finditer
from Bio import SeqIO
from Bio.Seq import Seq

usage = "Usage: python gap_finder.py genome.fa > n_windows.bed\n"
if len(sys.argv) < 2:
    sys.stderr.write("Need to specify a genome file.\n{}".format(usage))
    sys.exit(1)
if len(sys.argv) > 2:
    sys.stderr.write("Too many arguments.\n{}".format(usage))
    sys.exit(1)
if not os.path.isfile(sys.argv[1]):
    sys.stderr.write("File {} not found.\n{}".format(sys.argv[1], usage))
    sys.exit(1)

for record in SeqIO.parse(sys.argv[1], "fasta"):
    s = str(record.seq)

    for gap in finditer("N+", s):
        print(record.id, gap.start(), gap.end(), sep="\t")
