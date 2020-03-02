import sys
import os
from re import finditer
from Bio import SeqIO

usage = "Usage: python gap_finder.py genome.fa > n_windows.bed\n"
if len(sys.argv) < 2:
    sys.stderr.write("Need to specify a genome file.\n{}".format(usage))
    sys.exit(1)
if len(sys.argv) > 3:
    sys.stderr.write("Too many arguments.\n{}".format(usage))
    sys.exit(1)
if not os.path.isfile(sys.argv[1]):
    sys.stderr.write("File {} not found.\n{}".format(sys.argv[1], usage))
    sys.exit(1)

try:
    min_gap_length = sys.argv[2]
except IndexError:
    min_gap_length = str(1)

# Set count_only = False to print gap ranges to standard out
# and messages to standard error (default behavior)

# Set count_only = True to report number of gaps found but omit ranges
# (useful to test different minimum gap sizes).
count_only = False

sys.stderr.write("Finding gaps of minimum length {} in file {}\n".format(min_gap_length, sys.argv[1]))
if count_only: sys.stderr.write("Not printing intervals; counting gaps only\n")

gap_count = 0

for record in SeqIO.parse(sys.argv[1], "fasta"):
    s = str(record.seq)
    sys.stderr.write("Counting gaps in record {}\n".format(record.id))

    for gap in finditer("N{{{},}}".format(min_gap_length), s):
        if not count_only:
            print(record.id, gap.start(), gap.end(), sep="\t")
        gap_count += 1

sys.stderr.write("Total gaps found: {}\n".format(gap_count))
