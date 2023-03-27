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
parser.add_option("--names", dest="NA", help="Output file")
parser.add_option("--output", dest="OUT", help="Output file")
parser.add_option("--pools", dest="POO", help="Output file")
parser.add_option("--exclude", dest="EX", help="logical parameter")

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


Poslist = [0, 1, 2]
EX = d(str)
Names = []
Poo = []

for l in load_data(options.EX):
    EX[l.rstrip()]

C = 3
OUT1 = open(options.OUT+"_names_sub.txt", "wt")
for l in load_data(options.NA):
    if l.rstrip() in EX:
        C += 1
        continue
    Poslist.append(C)
    C += 1
    OUT1.write(l)

OUT2 = open(options.OUT+"_pools_sub.txt", "wt")
C = 3
for l in load_data(options.POO):
    if C not in Poslist:
        C += 1
        continue
    C += 1
    OUT2.write(l)

for l in load_data(options.IN):
    a = l.rstrip().split()
    print("\t".join([a[i] for i in Poslist]))

OUT1.close()
OUT2.close()
