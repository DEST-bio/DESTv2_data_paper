import sys
from collections import defaultdict as d
from optparse import OptionParser, OptionGroup

#Author: Martin Kapun

#########################################################   HELP   #########################################################################
usage="python %prog --input file --output file "
parser = OptionParser(usage=usage)
group=OptionGroup(parser,'< put description here >')

#########################################################   CODE   #########################################################################

parser.add_option("--input", dest="IN", help="Input file")
parser.add_option("--MAF", dest="MAF", help="GLOBAL minimum allele frequency thresholds as comma separated list")
parser.add_option("--ID", dest="ID", help="optionally Sample ID's",default="na")

(options, args) = parser.parse_args()
parser.add_option_group(group)

def load_data(x):
  ''' import data either from a gzipped or or uncrompessed file or from STDIN'''
  import gzip
  if x=="-":
      y=sys.stdin
  elif x.endswith(".gz"):
      y=gzip.open(x,"rt", encoding="latin-1")
  else:
      y=open(x,"r", encoding="latin-1")
  return y

NS,SS=d(lambda: d(int)),d(lambda: d(int))
stepsize=[float(x) for x in options.MAF.split(",")]
C=1
for l in load_data(options.IN):
    #print(l)
    if l.startswith("##"):
        continue
    a=l.rstrip().split()
    if l.startswith("#"):
        header=a[9:]
        if options.ID!="na":
            Keep=[x for x in range(len(header)) if header[x] in options.ID.split(",") ]
        continue
    INFO=a[7]
    RC=float(INFO.split("=")[1].split(";")[0])
    AC=float(INFO.split("=")[2].split(";")[0])
    if RC+AC==0:
        continue
    if AC/(AC+RC)<=0.5:
        MA=AC/(AC+RC)
    else:
        MA=1-AC/(AC+RC)
    #print(RC,AC,MA,round(MA/float(options.SS),0))
    if options.ID!="na":
        GT=[a[9:][x] for x in Keep if a[9:][x] != "./."]

        if len(list(set(GT)))==1:
            #print(GT)
            continue
    for x in stepsize:
        if MA<x:
            continue
        if "missense_variant" in INFO:
            NS[x][a[0]]+=1
        if "synonymous_variant"in INFO:
            SS[x][a[0]]+=1

for i in sorted(NS.keys()):
    ns,ss=0,0
    for k in sorted(NS[i].keys()):
        print(str(i),k,str(NS[i][k]),str(SS[i][k]),str(NS[i][k]/SS[i][k]),sep='\t')
        ns+=NS[i][k]
        ss+=SS[i][k]
    print(str(i),"genomewide",str(ns),str(ss),str(ns/ss),sep='\t')
