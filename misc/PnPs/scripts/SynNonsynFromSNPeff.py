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


for l in load_data(options.IN):
    Type = []
    a = l.rstrip().split()
    if l.startswith("#"):
        continue
    if "synonymous_variant" in a[7]:
        Type.append("Syn")
    elif "missense_variant" in a[7]:
        Type.append("Non")
    if Type == [] or len(Type) == 2:
        continue
    print("\t".join([a[0], a[1], ",".join(Type)]))
