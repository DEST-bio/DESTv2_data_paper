# Authors: M Kapun & S Steindl
# Analysis of Inversion Markers of the DESTv2 data


# data availability: http://berglandlab.uvadcos.io/gds/ 
# https://github.com/DEST-bio/DESTv2/blob/main/populationInfo/dest_v2.samps_25Feb2023.csv



# Convert VCF to SYNC

## PoolSNP
gunzip -c /media/DEST2_NHM/data/dest.all.PoolSNP.001.50.25Feb2023.norep.vcf.gz |
    parallel \
        --jobs 80 \
        --pipe \
        -k \
        --cat python3 /media/DEST2_NHM/scripts/VCF2sync.py \
        --input {} |
    gzip >/media/DEST2_NHM/output/dest.all.PoolSNP.001.50.25Feb2023.norep.gz

## additional command repetition for SNAPE 
gunzip -c /media/DEST2_NHM/data/dest.PoolSeq.SNAPE.NA.NA.25Feb2023.norep.vcf.gz |
    parallel \
        --jobs 80 \
        --pipe \
        -k \
        --cat python3 /media/DEST2_NHM/scripts/VCF2sync.py \
        --input {} |
    gzip >/media/DEST2_NHM/output/dest.PoolSeq.SNAPE.NA.NA.25Feb2023.norep.gz


# Get count at inversion specific marker SNPs (repeat for eather PoolSNP or SNAPE)
gunzip -c /media/DEST2_NHM/output/dest.PoolSeq.SNAPE.NA.NA.25Feb2023.norep.gz |
    parallel \
        --pipe \
        --jobs 20 \
        -k \
        --cat python3 /media/DEST2_NHM/scripts/overlap_in_SNPs.py \
        --source /media/DEST2_NHM/scripts/inversion_markers_v6.txt_pos \
        --target {} \
        >output/SNAPE_inversion_markers.sync

### check if variable call works
NAMES=$(gunzip -c /media/DEST2_NHM/data/dest.PoolSeq.SNAPE.NA.NA.25Feb2023.norep.vcf.gz | head -50 | awk '/^#C/' | cut -f10- | tr '\t' ',')

# Calculate average frequencies for marker SNPs
python3 /media/DEST2_NHM/scripts/inversion-freqs.py \
    /media/DEST2_NHM/scripts/inversion_markers_v6.txt \
    /media/DEST2_NHM/output/SNAPE_inversion_markers.sync \
    $NAMES \
    >/media/DEST2_NHM/results/SNAPE_nhm_inversion.af

# Plot Correlation of Inversion Markers of PoolSNP & SNAPE 
# Rscript /media/DEST2_NHM/scripts/PlotCorr.R


# Plot for each Inversion occuring in the data set
# RScript "LatContInv.R" uses three input arguments (AF table, metadata, name of inversion)
inversions=("In.2L.t" "In.2R.Ns" "In.3L.P" "In.3R.C" "In.3R.K" "In.3R.Mo" "In.3R.Payne")

for i in ${inversions[@]}; do
    inv=$(expr $i)
    #echo $inv
    Rscript /media/DEST2_NHM/scripts/LatContInv.R /media/DEST2_NHM/results/PoolSNP_nhm_inversion.af /media/DEST2_NHM/data/dest_v2.samps_25Feb2023.csv ${inv}
done



##Additional script to directly work in R and make adjustments:
## /media/DEST2_NHM/scripts/PlotInversions.R

## Get Coverages (DP) from sync files as weights for regression analysis
#python3 /media/DEST2_NHM/scripts/SubsampleSyncCov.py --sync PoolSNP_inversion_markers.sync  > CoveragesPoolSNP.csv
python3 /media/DEST2_NHM/scripts/SubsampleSyncCov.py --sync SNAPE_inversion_markers.sync  > CoveragesSNAPE.csv


python3 /media/DEST2_NHM/scripts/markerfreqs.py \
    /media/DEST2_NHM/scripts/inversion_markers_v6.txt \
    /media/DEST2_NHM/output/SNAPE_inversion_markers.sync \
    $NAMES > /media/DEST2_NHM/output/SNAPE_singlefreqs.txt

python3 /media/DEST2_NHM/scripts/markerfreqs.py \
/media/DEST2_NHM/scripts/inversion_markers_v6.txt \
/media/DEST2_NHM/output/PoolSNP_inversion_markers.sync \
$NAMES > /media/DEST2_NHM/output/PoolSNP_singlefreqs.txt

python3 /media/DEST2_NHM/scripts/CovPerMarker.py --sync SNAPE_inversion_markers.sync  > SNAPE_CovPerMarker.csv
