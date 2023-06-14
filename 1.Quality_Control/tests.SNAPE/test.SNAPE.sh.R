#### 
#### 

SNAPE=/scratch/yey2sn/SNAPE_tinkering/snape-pooled/snape-pooled

mpile=/project/berglandlab/DEST/dest_mapped/RECENT_OUTPUTS/US_Tex_Lub_0_2013-09-15/US_Tex_Lub_0_2013-09-15.mel_mpileup.txt.gz

zcat $mpile | head -n 25000 > test.mpileup

####
####
awk '{if (last != $1) close(last); print >> $1; last = $1}' test.mpileup

nflies=100
chr=2L
theta=0.005
D=0.01
priortype="informative"
fold="unfolded"
sample="test"

$SNAPE -nchr $(($nflies*2)) -theta $theta -D $D -priortype $priortype -fold $fold < ${chr} > ${chr}-$sample-SNAPE.txt
