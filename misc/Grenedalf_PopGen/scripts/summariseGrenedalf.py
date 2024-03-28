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
parser.add_option("--coverage", dest="COV", help="Output file")
parser.add_option("--meta", dest="MET",
                  help="logical parameter")
parser.add_option("--chrom", dest="chromosome",
                  help="logical parameter", action="store_true")


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


Chrom = ("4", "2L", "2R", "3L", "3R", "X")
Lengths = (1348131, 23513712, 25286936, 28110227, 32079331, 23542271)

CHROMDICT = {k: v for k, v in zip(*[Chrom, Lengths])}

if options.chromosome:
    Cov = d(lambda: d(str))

    for l in load_data(options.COV):
        a = l.rstrip().split()
        if a[1] not in CHROMDICT:
            continue
        Cov[a[0]][a[1]] = CHROMDICT[a[1]]-int(a[4])
else:
    Cov = d(lambda: d(lambda: d(str)))

    for l in load_data(options.COV):
        a = l.rstrip().split()
        Cov[a[0]][a[1]][(int(a[2])+int(a[3]))/2] = int(a[3]) - \
            int(a[2])-1-int(a[4])

header = ""
META = d(lambda: d())
for l in load_data(options.MET):
    a = l.rstrip().split(",")
    if header == "":
        header = a
        continue
    DATA = {k: v for k, v in zip(*[header, a])}
    META[a[0]]["continent"] = DATA["continent"]
    META[a[0]]["country"] = DATA["country"]
    META[a[0]]["lat"] = DATA["lat"]
    META[a[0]]["long"] = DATA["long"]

header = ""
for l in load_data(options.IN):
    a = l.rstrip().split(",")
    chrom, start, end = a[:3]
    pops = a[3:]
    if header == "":
        header = pops
        continue
    DATA = {k: v for k, v in zip(*[header, pops])}
    for ID, v in sorted(DATA.items()):
        # print(ID)
        k, Type = ID.split(".")

        if k not in META:
            continue
        if k not in Cov:
            continue
        if chrom not in Cov[k]:
            continue
        # get denominator for windowwise averages
        if options.chromosome:
            REL = str(float(v)/Cov[k][chrom])
        else:
            if (int(start)+int(end))/2 not in Cov[k][chrom]:
                REL = "NA"
            elif Cov[k][chrom][(int(start)+int(end))/2] == 0:
                REL = "NA"
            else:
                REL = str(float(v)/Cov[k][chrom][(int(start)+int(end))/2])
        if "abs" in Type:
            print(k,
                  META[k]["continent"],
                  META[k]["country"],
                  META[k]["lat"],
                  META[k]["long"],
                  chrom,
                  str((int(start)+int(end))/2),
                  Type.replace("abs", "rel"),
                  REL,
                  sep="\t")
        if "tajimas_d" in Type:
            print(k,
                  META[k]["continent"],
                  META[k]["country"],
                  META[k]["lat"],
                  META[k]["long"],
                  chrom,
                  str((int(start)+int(end))/2),
                  Type.replace("abs", "rel"),
                  v,
                  sep="\t")
