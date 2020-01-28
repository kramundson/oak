# Lisa Malins
# 24 Jan 2020
# GC_content_summarizer.py
"""
Summarizes frequency of G, C, and N within windows of a genome.
Usage: python gcn_content_summarizer.py windows.bed genome.fa > output.tsv
"""

import sys
from Bio import SeqIO

# Create dictionary of window regions
def parse_windows(windows_file):
    windows = {}

    with open(windows_file, 'r') as f:
        while True:
            line = f.readline()
            if not line: break

            l = line.split()
            name = l[0]

            # Create dictionary entry for sequence if does not exist
            if name not in windows:
                windows[name] = []

            # Save region into dictionary
            windows[name].append([int(i) for i in l[1:3]])

    return windows

def main():

    # Parse file of windows into dictionary
    windows_file = sys.argv[1]
    windows = parse_windows(windows_file)

    # Loop through sequences in genome fasta,
    # count G/C/N content in each region,
    # and write out in tab-delimited format.

    genome = sys.argv[2]
    with open("/dev/fd/1", 'w') as out:

        # Loop through records in fasta
        for record in SeqIO.parse(genome, "fasta"):

            # Loop through regions
            for region in windows[record.id]:
                # Save start and end as integers
                start = region[0]
                end = region[1]

                # Save counts of G, C, and N
                counts = [str(record.seq.count(l, start, end)) for l in ["G", "C", "N"]]

                # Print sequence name, start, end, and G/C/N counts
                print("\t".join([record.id, str(start), str(end)] + counts))



if __name__ == "__main__":
    main()
