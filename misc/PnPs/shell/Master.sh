#!/bin/sh

#  Master.sh
#
#
#  Created by Martin Kapun on 04.11.20.
#

### define syn/non-syn SNPs

mkdir /Volumes/MartinResearch2/DEST/analyses/PnPs

gunzip -c /Volumes/MartinResearch2/DEST/FinalData/dest.PoolSeq.SNAPE.001.50.ann.vcf.gz \
| parallel \
-k \
--pipe \
-j 22 \
--no-notice \
--cat python /Users/mkapun/Documents/GitHub/UsefullStuff/UsefulScripts/SynNonsynFromSNPeff.py \
--input {} | gzip > /Volumes/MartinResearch2/DEST/analyses/PnPs/snape_snps.pnps.gz


gunzip -c /Volumes/MartinResearch2/DEST/FinalData/dest.July6_2020.001.10.ann.vcf.gz \
| parallel \
-k \
--pipe \
-j 22 \
--no-notice \
--cat python /Users/mkapun/Documents/GitHub/UsefullStuff/UsefulScripts/SynNonsynFromSNPeff.py \
--input {} | gzip > /Volumes/MartinResearch2/DEST/analyses/PnPs/poolsnp_snps.pnps.gz


python3 /Users/mkapun/Documents/GitHub/DEST/utils/PNPS4DGRP.py \
--input /Volumes/MartinResearch1/TadExpEvol/new/data/dgrp2_ann.vcf.gz \
--MAF 0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2 \
> /Volumes/MartinResearch2/DEST/analyses/PnPs/null.pnps

python3 /Users/mkapun/Documents/GitHub/DEST/utils/PNPS4SNAPE.py \
--input /Volumes/MartinResearch2/DEST/data/dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz \
--output /Volumes/MartinResearch2/DEST/analyses/PnPs/SNAPE_full.pnps.gz \
--MAF 0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2

python3 /Users/mkapun/Documents/GitHub/DEST/utils/PNPS4SNAPE.py \
--input /Volumes/MartinResearch2/DEST/data/dest.PoolSeq.SNAPE.NA.NA.10Nov2020.ann.vcf.gz \
--output /Volumes/MartinResearch2/DEST/analyses/PnPs/SNAPE_full.pnps.gz \
--MAF 0


### now also for parameters:

mkdir /Volumes/MartinResearch2/DEST/analyses/PnPs/ParamPoolSNP

cd /Volumes/MartinResearch2/DEST/data/paramTest

for i in *.gz

do

  python3 /Users/mkapun/Documents/GitHub/DEST/utils/PNPS4SNAPE.py \
  --input ${i} \
  --output /Volumes/MartinResearch2/DEST/analyses/PnPs/ParamPoolSNP/${i}.txt \
  --Parse \
  --MAF 0.0 &

done

printf "MAF\tMAC\tMAFc\tPOP\tChrom\tNS\tSS\tpNpS\n" > /Volumes/MartinResearch2/DEST/analyses/PnPs/ParamPoolSNP.txt

for i in /Volumes/MartinResearch2/DEST/analyses/PnPs/ParamPoolSNP/*.gz.txt

do

gunzip -c $i | awk 'NR>1' >> /Volumes/MartinResearch2/DEST/analyses/PnPs/ParamPoolSNP.txt

done
