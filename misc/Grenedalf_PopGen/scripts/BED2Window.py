import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import math
# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(
    parser, 'This script calculates the number of sites per window that do not pass the quality criteria using the BED output of the DEST pipeline script')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--window", dest="W", help="Output file")

(options, args) = parser.parse_args()
parser.add_option_group(group)


def load_data(x):
    ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
    import gzip
    if x == "-":
        y = sys.stdin
    elif x.endswith(".gz"):
        y = gzip.open(x, "rt", encoding="latin-1")
    else:
        y = open(x, "r", encoding="latin-1")
    return y


INPUT = options.IN

if options.W == "chrom":
    WINDOW = "chrom"
else:
    WINDOW = int(options.W)

NAME = INPUT.split("/")[-1].split(".bed")[0]

if WINDOW != "chrom":
    DEL = d(lambda: d(int))

    for l in load_data(INPUT):
        C, S, E = l.rstrip().split()
        for i in range(int(S)+1, (int(E)+1)):
            DEL[C][int(math.floor(i/WINDOW)*WINDOW)] += 1

    for chrom, v in sorted(DEL.items()):
        for k, v1 in sorted(v.items()):
            print(NAME, chrom, str(k+1), str(k+WINDOW), str(v1), str(v1/WINDOW))

else:
    DEL = d(int)

    for l in load_data(INPUT):
        C, S, E = l.rstrip().split()
        for i in range(int(S)+1, (int(E)+1)):
            DEL[C] += 1

    for chrom, v in sorted(DEL.items()):
        print(NAME, chrom, "NA\tNA", str(v), str(v))
