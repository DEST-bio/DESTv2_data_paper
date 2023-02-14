### prepare hologenome
cat \
holo_dmel_6.12.fa \
D_simulans.fasta > \
dmel_holo_sim.fa

bwa index dmel_holo_sim.fa


module load gcc/9.2.0
module load bwa/0.7.17

bwa index dmel_holo_sim.fa