# ijob -A berglandlab_standard -c1 -p dev --mem=1G
# module load gcc/7.1.0 openmpi/3.1.4 R/3.6.3; R

library(SeqArray)

new_geno <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds")
old_geno <- seqOpen("/project/berglandlab/DEST/gds/dest.all.PoolSNP.001.50.26April2023.norep.ann.gds")

seqSetFilter(new_geno, sample.id=c("ES_Bar_Bar_1_2012-07-13", "ES_Bar_Bar_1_2012-10-16",
                                   "AT_Wie_Gro_1_2012-08-03", "AT_Wie_Gro_1_2012-10-20"), variant.id=10)

seqSetFilter(old_geno, sample.id=c("ES_Bar_Bar_1_2012-07-13", "ES_Bar_Bar_1_2012-10-16",
                                  "AT_Wie_Gro_1_2012-08-03", "AT_Wie_Gro_1_2012-10-20"), variant.id=10)


seqGetData(new_geno,"annotation/format/DP")
seqGetData(old_geno,"annotation/format/DP")
