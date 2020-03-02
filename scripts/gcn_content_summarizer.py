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
    sys.stderr.write("Parsing regions from {}\n".format(windows_file))

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

    sys.stderr.write("Done parsing regions from {}\n\n".format(windows_file))

    return windows

def main():
    usage="Usage: python gcn_content_summarizer.py windows.bed genome.fa > output.tsv\n"
    if len(sys.argv) != 3:
        sys.stderr.write("Wrong number of arguments.\n" + usage)
        sys.exit(1)

    sys.stderr.write("Windows file: {} \t\tGenome file: {}\n\n".format(sys.argv[1], sys.argv[2]))

    # Parse file of windows into dictionary
    windows_file = sys.argv[1]
    windows = parse_windows(windows_file)

    # Loop through sequences in genome fasta,
    # count G/C/N content in each region,
    # and write out in tab-delimited format.
    sys.stderr.write("Counting G/C/N content in genome {}\n".format(sys.argv[2]))

    genome = sys.argv[2]
    with open("/dev/fd/1", 'w') as out:
        # Write column names for R
        print("\t".join(["chrom", "chromStart", "chromEnd", "G", "C", "N"]))

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

    sys.stderr.write("Done counting G/C/N content in genome {}\n".format(sys.argv[2]))


if __name__ == "__main__":
    main()
