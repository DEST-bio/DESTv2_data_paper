import pandas as pd
import gzip

"""
This script is to be run before get_average_n_eff_per_samp.py, which uses the output
of this script to calculate per-popinfo average effective sample sizes (n_eff) needed
for allele count rounding during subsequent site frequency spectrum construction
with `construct_sfs.py`.
"""

# This script is to be run before get_average_n_eff_per_popinfo.py

def main():
    # Setup
    data_dir = "/data/rmccoy22/dbass13/DEST2/"
    vcf_file = data_dir + "vcfs/dest2_clustered_mainChroms.vcf.gz"
    metadata_file = data_dir + "metadata/dest_v2.samps_8Jun2023.csv"
    output_file = data_dir + "metadata/mean_n_eff_per_sample_on_mainChroms.csv"

    get_average_n_eff_per_samp(vcf_file, metadata_file, output_file)

def get_average_n_eff_per_samp(vcf_file: str, 
                               metadata_file: str, 
                               output_file: str) -> None:
    """
    :param vcf_file str: Name of VCF file, including path
    :param metadata_file str: Name of metadata file, including path. The file is sourced
        from https://github.com/DEST-bio/DESTv2/blob/main/populationInfo/dest_v2.samps_8Jun2023.csv
    :param output_file str: Name of output file, including path
    :return: None, but writes to an output file the name of the sample and the
        corresponding mean effective coverage of the sample
    :rtype: None
    """
    
    # Load in DEST 2.0 metadata into samp_dict, which maps each sampleID to
    # a 3-array of [nFlies, sum(n_eff), num_sites]
    metadata = pd.read_csv(metadata_file)
    samp_to_nSamp = pd.Series(metadata.nFlies.values,
                              index=metadata.sampleId).to_dict()
    samp_dict = {samp: [nFlies, 0, 0] 
                 for samp, nFlies in samp_to_nSamp.items()}
    
    # Load VCF file
    vcf = gzip.open(vcf_file)
    depth_info_index = 3 # index of DP in the FORMAT field
    
    # Iterate over lines of VCF file
    for line in vcf:
        line = line.decode() # Account for any UTF-8 wackiness

        # Skip intro info
        if line.startswith("##"):
            continue
        # Read sample IDs from header
        if line.startswith("#"):
            sampleIds = line.split()[9:]
            continue
        
        # Iterate over samples in line
        for i, field in enumerate(line.split()[9:]):
            # Skip missing sequencing data
            if field.startswith("."):
                continue
            nFlies = samp_dict[sampleIds[i]][0]
            nChromosomes = nFlies * 2 # diploid
            depth = int(field.split(":")[depth_info_index])
            samp_dict[sampleIds[i]][1] += n_eff(nChromosomes, depth)
            samp_dict[sampleIds[i]][2] += 1

    vcf.close()

    # Take averages of n_eff over all sites included in each sample's column.
    # Many samples will have a mean n_eff of -1, which indicates that no sites
    # were included in the sample's column, but these are filtered out as a
    # convenient side effect of selecting for the sample within each locality
    # with greatest nFlies in the R script that generates popinfos.
    samp_to_mean_n_eff = {samp: arr[1] / arr[2] if arr[2] > 0 else -1
                          for samp, arr in samp_dict.items()}

    # Write output as CSV
    samp_to_mean_n_eff_df = pd.DataFrame(samp_to_mean_n_eff.items(), 
                                         columns=['sampleId', 'mean_n_eff'])
    samp_to_mean_n_eff_df.to_csv(output_file, index=False)

def n_eff(nChromosomes: int, 
          depth: int) -> float:
    """
    :param nChromosomes int: Number of chromosomes sampled from the population,
        equals 2 x n for diploid individuals
    :param depth int: Average sequencing coverage across the region-of-interest
        of the genome
    :return: Effective sample size (n_eff) for the sample, which tunes down the
        sequencing depth in order to account for uncertainty introduced by pool-seq
    :rtype: float
    """

    return (nChromosomes * depth - 1) / (nChromosomes + depth)

if __name__ == "__main__":
    main()