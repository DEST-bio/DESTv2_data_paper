import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup
import gzip

# Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage = "python %prog --input file --output file "
parser = OptionParser(usage=usage)
group = OptionGroup(parser, '< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--output", dest="OUT", help="Output file")
parser.add_option("--offset", dest="OFF", help="offset to breakpoints in bp")
parser.add_option("--logical", dest="log",
                  help="logical parameter", action="store_true")
parser.add_option("--param", dest="param",
                  help="numerical parameter", default=1)

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


Inv = d(list)
InvOut = d(str)

for l in load_data(options.param):
    a = l.rstrip().split(",")
    for i in range(int(a[2])-int(options.OFF), int(a[3])+int(options.OFF)):
        Inv[a[1]+str(i)].append(a[0])
        # print(a[1]+str(i))
    InvOut[a[0]] = gzip.open(options.OUT+"_"+a[0]+".af.gz", "wt")

Out = gzip.open(options.OUT+"_NoInv.af.gz", "wt")

for l in load_data(options.IN):
    a = l.rstrip().split()
    if a[0]+a[1] in Inv:
        for inv in Inv[a[0]+a[1]]:
            InvOut[inv].write(l)
    else:
        Out.write(l)
