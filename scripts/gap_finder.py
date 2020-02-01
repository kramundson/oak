import sys
from re import finditer
from Bio import SeqIO
from Bio.Seq import Seq

usage = "python gap_finder.py genome.fa n_windows.bed"

for record in SeqIO.parse(sys.argv[1], "fasta"):
    s = str(record.seq)

    for gap in finditer("N+", s):
        print(record.id, gap.start(), gap.end(), sep="\t")
