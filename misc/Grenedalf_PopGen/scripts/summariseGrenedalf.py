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
parser.add_option("--meta", dest="MET",
                  help="logical parameter")


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


FULL = d(lambda: d(lambda: d(float)))
print("ID\tContinent\tCountry\tChrom\tStat\tValue")
for i in range(len(Chrom)):
    header = ""
    for l in load_data(options.IN+"_"+Chrom[i]+"diversity.csv"):
        a = l.rstrip().split()
        if header == "":
            header = a
            continue
        DATA = {k: v for k, v in zip(*[header, a])}
        for k, v in sorted(META.items()):
            for Type in ["theta_pi_abs", "theta_watterson_abs", "tajimas_d", "snp_count"]:
                if k+"."+Type not in DATA:
                    continue

                if "abs" in Type:
                    if DATA[k+"."+Type] == "NA":
                        Value = "NA"
                    else:
                        Value = str(float(DATA[k+"."+Type])/Lengths[i])
                else:
                    Value = DATA[k+"."+Type]
                print(k, META[k]["continent"], META[k]
                      ["country"], Chrom[i], Type, Value, sep="\t")
                if DATA[k+"."+Type] != "NA":
                    FULL[k][Type]["counts"] += float(Value)*Lengths[i]
                    FULL[k][Type]["lengths"] += Lengths[i]

for k, v in sorted(FULL.items()):
    for Type in ["theta_pi_abs", "theta_watterson_abs", "tajimas_d", "snp_count"]:
        if v[Type]["lengths"] == 0:
            print
            WAv = "NA"
        else:
            WAv = str(v[Type]["counts"]/v[Type]["lengths"])
        print(k, META[k]["continent"], META[k]
              ["country"], "GenomeWide", Type, WAv, sep="\t")
