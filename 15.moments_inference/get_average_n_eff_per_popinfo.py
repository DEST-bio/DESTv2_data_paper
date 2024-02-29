import pandas as pd
import sys

"""
This script is to be run after get_average_n_eff_per_samp.py in order to calculate
average effective sample sizes (n_eff) for the samples listed in a single popinfo
file before constructing site frequency spectra from those popinfos with `construct_sfs.py`.
"""

def main():
    # Setup
    data_dir = "/project/berglandlab/david_DEST2/"
    metadata_file = data_dir + "metadata/mean_n_eff_per_sample_on_mainChroms.csv"
    popinfo_name = sys.argv[1]
    output_file = data_dir + "metadata/mean_n_eff_per_popinfo.csv"

    get_average_n_eff_per_popinfo(data_dir, metadata_file, popinfo_name, output_file)

def get_average_n_eff_per_popinfo(data_dir: str, 
                                  metadata_file: str, 
                                  popinfo_name: str, 
                                  output_file: str) -> None:
    """
    :param data_dir str: Path to directory containing data, a directory `popinfos/`
        that contains popinfo files 
    :param metadata_file str: Name of metadata, including path. The file is sourced
        from https://github.com/DEST-bio/DESTv2/blob/main/populationInfo/dest_v2.samps_8Jun2023.csv
    :param popinfo_name str: Name of input popinfo file, no path
    :param output_file str: Name of output file, including path
    :return: None, but appends to an output file the name of the popinfo file and
        the corresponding mean effective coverage of the samples listed in the popinfo
        file
    :rtype: None
    """

    # Load in data
    popinfo_file = data_dir + "popinfos/" + popinfo_name
    metadata_df = pd.read_csv(metadata_file)
    popinfo = pd.read_csv(popinfo_file, sep="\t", header=None,
                          names=["sampleId", "popId"])

    # Calculate average effective coverage
    mean_n_eff = popinfo.merge(metadata_df, 
                               on='sampleId', 
                               how='left').mean_n_eff.mean()
    print(mean_n_eff)
    
    # Write output
    with open(output_file, 'a') as f:
        f.write(popinfo_name + "," + str(mean_n_eff) + "\n")

if __name__ == "__main__":
    main()
