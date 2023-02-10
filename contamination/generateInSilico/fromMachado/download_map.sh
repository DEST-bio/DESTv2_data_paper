#!/bin/bash
######## Install seqTk
#=======#   git clone https://github.com/lh3/seqtk.git;
#=======#   cd seqtk; make

######## Get the RunInfo tables for D. simulans and D. melanogaster reads for mapping simulation
    #### simalans SRA data
    ### https://www.ncbi.nlm.nih.gov/Traces/study/?WebEnv=NCID_1_11549260_130.14.22.76_5555_1552330190_2321860912_0MetA0_S_HStore&query_key=10
    ### RunInfo Table gets dumped to file : /mnt/pricey_4/sim_mel_contamination_simulation/metadata/simulans_SRA


    ### melanogaster SRA data
    ### https://www.ncbi.nlm.nih.gov/Traces/study/?acc=DGRP&go=go
    ### table gets dumped to file : /mnt/pricey_4/sim_mel_contamination_simulation/metadata/melanogaster_SRA


######## Specify which SRA files to download
    howManySRA=1

    rm /mnt/pricey_4/sim_mel_contamination_simulation/metadata/toDownLoad.delim
    grep "HiSeq" < /mnt/pricey_4/sim_mel_contamination_simulation/metadata/simulans_SRA | cut -f3,6 | head -n${howManySRA} | sed 's/$/\tsim/g' > /mnt/pricey_4/sim_mel_contamination_simulation/metadata/toDownLoad.delim
    grep "HiSeq" < /mnt/pricey_4/sim_mel_contamination_simulation/metadata/melanogaster_SRA | cut -f16,19 |awk '{print $2"\t"$1"\t"$3}' | head -n${howManySRA} | sed 's/$/\tmel/g' >> /mnt/pricey_4/sim_mel_contamination_simulation/metadata/toDownLoad.delim


######## Download raw sequence data
    downloadSRA () {
            echo ${1}

            ### pull out SRR number
            srr=$( echo ${1} | cut -f2 -d' ' )
            fileStem=$( echo ${1} | cut -f1,3 -d' ' | tr ' ' '_' )


            ### srr=SRR3585384
            ### fileStem=foo
            firstThree=${srr:0:6}

            #pop=$( echo $outDir | cut -f1 -d'/' )

            #if [ ! -d /mnt/inbred/inbredLines/00_rawData/${pop} ]; then
            #       mkdir /mnt/inbred/inbredLines/00_rawData/${pop}
            #fi
            #
            #if [ ! -d /mnt/inbred/inbredLines/00_rawData/${outDir} ]; then
            #       mkdir /mnt/inbred/inbredLines/00_rawData/${outDir}
            #fi

           wget ftp://ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${firstThree}/${srr}/${srr}.sra \
           -O /mnt/pricey_4/sim_mel_contamination_simulation/sra/${fileStem}.sra

            #/home/bergland/.aspera/connect/bin/ascp -T -k 1 \
            #-i /home/bergland/.aspera/connect/etc/asperaweb_id_dsa.openssh \
            #anonftp@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${firstThree}/${srr}/${srr}.sra \
            #/mnt/pricey_4/sim_mel_contamination_simulation/sra/${fileStem}.sra

            fastq-dump --split-files \
            --gzip \
            -O /mnt/pricey_4/sim_mel_contamination_simulation/fastq \
            /mnt/pricey_4/sim_mel_contamination_simulation/sra/${fileStem}.sra

            #rm /mnt/inbred/inbredLines/00_rawData/${outDir}/${srr}.sra

    }
    export -f downloadSRA

   # nohup parallel --gnu -j1 -a /mnt/pricey_4/sim_mel_contamination_simulation/metadata/toDownLoad.delim downloadSRA &

######## Downsample reads
    downSample_mix () {

        fracReads_sim=${1}    #; #fracReads_sim=.1
        fracReads_mel=$( echo 1-${fracReads_sim} | bc -l )

        nTargetReads=${2} #; #nTargetReads=10001

        nReads_sim=$( echo  "(${nTargetReads}*${fracReads_sim}+0.5)/1" | bc )
        nReads_mel=$( echo  "(${nTargetReads}*${fracReads_mel}+0.5)/1" | bc )

        echo ${nReads_sim}

        /mnt/pricey_4/sim_mel_contamination_simulation/software/seqtk/seqtk sample \
        -s seed=11 -2 \
        /mnt/pricey_4/sim_mel_contamination_simulation/fastq/Sz100_sim_1.fastq.gz \
        ${nReads_sim} > /mnt/pricey_4/sim_mel_contamination_simulation/downSampleFastq/simDown_${nReads_sim}.fa

        /mnt/pricey_4/sim_mel_contamination_simulation/software/seqtk/seqtk sample \
        -s seed=11 -2 \
        /mnt/pricey_4/sim_mel_contamination_simulation/fastq/DGRP-395_mel_1.fastq.gz \
        ${nReads_mel} > /mnt/pricey_4/sim_mel_contamination_simulation/downSampleFastq/melDown_${nReads_mel}.fa



    }
    export -f downSample_mix

    #nohup parallel --gnu -j1 downSample_mix ::: 0 .1 .2 .3 .4 .5 .6 .7 .8 ::: 1000000 &
     parallel --gnu -j1 downSample_mix ::: .000001 .00001 .0001 .001 .01 .1 ::: 1000000


######## generate composite genome
    wget -O /mnt/pricey_4/sim_mel_contamination_simulation/ref/dsim-all-chromosome-r2.01.fasta.gz \
    ftp://ftp.flybase.net/genomes/Drosophila_simulans/dsim_r2.01_FB2015_01/fasta/dsim-all-chromosome-r2.01.fasta.gz

    wget -O /mnt/pricey_4/sim_mel_contamination_simulation/ref/dmel-all-chromosome-r5.5.fasta.gz \
    ftp://ftp.flybase.net//genomes/Drosophila_melanogaster/dmel_r5.5_FB2008_01/fasta/dmel-all-chromosome-r5.5.fasta.gz

    gunzip /mnt/pricey_4/sim_mel_contamination_simulation/ref/*

    sed -i 's/>/>dmel_/g' /mnt/pricey_4/sim_mel_contamination_simulation/ref/dmel-all-chromosome-r5.5.fasta
    sed -i 's/>/>dsim_/g' /mnt/pricey_4/sim_mel_contamination_simulation/ref/dsim-all-chromosome-r2.01.fasta

    cat \
    /mnt/pricey_4/sim_mel_contamination_simulation/ref/dmel-all-chromosome-r5.5.fasta \
    /mnt/pricey_4/sim_mel_contamination_simulation/ref/dsim-all-chromosome-r2.01.fasta > \
    /mnt/pricey_4/sim_mel_contamination_simulation/ref/combo.fasta

    bwa index /mnt/pricey_4/sim_mel_contamination_simulation/ref/combo.fasta

######## mix reads & map to composite genome
 mix_and_map () {

     fracReads_sim=${1}; #fracReads_sim=.00001
     fracReads_mel=$( echo 1-${fracReads_sim} | bc -l )

     nTargetReads=${2}; #nTargetReads=1000000

     nReads_sim=$( echo  "(${nTargetReads}*${fracReads_sim}+0.5)/1" | bc )
     nReads_mel=$( echo  "(${nTargetReads}*${fracReads_mel}+0.5)/1" | bc )


     cat /mnt/pricey_4/sim_mel_contamination_simulation/downSampleFastq/simDown_${nReads_sim}.fa \
     /mnt/pricey_4/sim_mel_contamination_simulation/downSampleFastq/melDown_${nReads_mel}.fa > \
     /mnt/pricey_4/sim_mel_contamination_simulation/downSampleFastq/mix.simDown_${nReads_sim}.melDown_${nReads_mel}.fa

     bwa mem -t10 \
     /mnt/pricey_4/sim_mel_contamination_simulation/ref/combo.fasta \
     /mnt/pricey_4/sim_mel_contamination_simulation/downSampleFastq/mix.simDown_${nReads_sim}.melDown_${nReads_mel}.fa | \
     samtools view -b - > /mnt/pricey_4/sim_mel_contamination_simulation/mapped_reads/mix.simDown_${nReads_sim}.melDown_${nReads_mel}.bam


 }
 export -f mix_and_map

 nohup parallel --gnu -j1 mix_and_map ::: 0 .1 .2 .3 .4 .5 .6 .7 .8 ::: 1000000 &

 nohup parallel --gnu -j1 mix_and_map ::: 0.000001 0.00001 0.0001 0.001 0.01 0.1::: 1000000 &

 ######## summarize
 summarize () {
     file=${1} #file=mix.simDown_0.melDown_1000000.bam
     samtools view ${file} | cut -f 1,3,4 | awk -v fn=${file} '{print fn"\t"$0}' > ${file}.delim
 }
 export -f summarize
 parallel --gnu -j1 summarize ::: $( ls /mnt/pricey_4/sim_mel_contamination_simulation/mapped_reads/*.bam)

### collapse files and extract out autosomes

cat /mnt/pricey_4/sim_mel_contamination_simulation/mapped_reads/*.delim | grep -E "2L|2R|3L|3R" | grep -v "Het" > \
/mnt/pricey_4/sim_mel_contamination_simulation/summary/output.delim
