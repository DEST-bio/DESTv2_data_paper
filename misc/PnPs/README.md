## _P_<sub>n</sub>_P_<sub>s</sub> analysis

Following the approach in our DESTv.1 paper, I now calculated _P_<sub>n</sub>_P_<sub>s</sub> ratios for the SNAPE dataset and counted the number of private SNPs. Based on a threshold for _P_<sub>n</sub>_P_<sub>s</sub> and privates SNPs (i.e., mean(stat)*SD(Stat)*1.96), I categorized populations with basedcalling based on SNAPE as "Keep" and "Exclude". The corresponding table can be found [here](results/classify_pops.txt) and the corresponding figure is below.

![Figure](results/classify.png)

To carry out the corresponding analysis with PoolSNP, I would need VCF datasets based on various MAC and MAF thresholds

