import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--output", dest="OUT", help="Output file")

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


C = 1
PRIVATE = d(int)
for l in load_data(options.IN):
    if l.startswith("##"):
        continue
    a = l.rstrip().split()
    if l.startswith("#"):
        header = a[9:]
        continue
    GT = [x.split(":")[0] for x in a[9:]]
    if GT.count("0/1") == 1:
        PRIVATE[header[GT.index("0/1")]] += 1
    if C % 1000000 == 0:
        print(C, "SNPs processed")
    C += 1

OUT = open(options.OUT, "wt")
OUT.write("POP\tN\n")
for k, v in sorted(PRIVATE.items()):
    OUT.write("\t".join([k, str(v)])+"\n")
