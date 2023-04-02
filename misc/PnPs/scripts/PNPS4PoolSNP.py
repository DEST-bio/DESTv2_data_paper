import sys
import gzip
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
parser.add_option("--MAC", dest="MAC",
                  help="GLOBAL minimum allele count thresholds as comma separated list")

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
stepsize = [int(x) for x in options.MAC.split(",")]
PRIVATE = d(lambda: d(lambda: d(lambda: d(int))))
for l in load_data(options.IN):
    if l.startswith("##"):
        continue
    a = l.rstrip().split()
    if l.startswith("#"):
        header = a[9:]
        continue
    if C % 1000000 == 0:
        print(C, "SNPs processed")
    C += 1
    if "missense_variant" not in a[7] and "synonymous_variant" not in a[7]:
        continue
    AC = sum([int(x)
             for x in a[7].split(";AC=")[1].split(";")[0].split(",")])
    GT = [x.split(":")[0] for x in a[9:]]
    GT_i = [header[x] for x in range(len(GT)) if GT[x] == "0/1"]
    for s in stepsize:
        if AC < s:
            continue
        for i in GT_i:
            if "missense_variant" in a[7]:
                PRIVATE[s][i][a[0]]["NS"] += 1
                PRIVATE[s][i]["genomewide"]["NS"] += 1
            if "synonymous_variant" in a[7]:
                PRIVATE[s][i][a[0]]["SS"] += 1
                PRIVATE[s][i]["genomewide"]["SS"] += 1

OUT = gzip.open(options.OUT, "wt")
OUT.write("MAC\tPOP\tChrom\tNS\tSS\tpNpS\n")
for MAC, v in sorted(PRIVATE.items()):
    for POP, v1 in sorted(v.items()):
        for Chrom, v2 in sorted(v1.items()):
            if Chrom == "Y":
                continue
            OUT.write(
                "\t".join([str(MAC), POP, Chrom, str(v2["NS"]), str(v2["SS"]), str(v2["NS"]/v2["SS"])])+"\n")
