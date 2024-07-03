import numpy as np
import moments
import sys
import os
import functools
import operator

"""
This script constructs site frequency spectra (SFS) from A VCF file representing
pool-seq data and a popinfo file. 
It is run in parallel array jobs with the wrapper `wrapper_construct_sfs.slurm`.
"""

# Global dictionary for mapping model names to (1) the mean effective sample size
# of samples included in the modeled populations and (2) an ordered list of the 
# names of the populations included in the model.
# TODO: Make this automatically import from `mean_n_eff_per_popinfo.csv` to circumvent
# manual addition of n_eff values to this dictionary.
popinfo_metadata = {"k4_Europe": [round(31.0888),
                                    ["Europe_west", "Europe_east"]],
                    "k8_Europe": [round(31.0888),
                                    ["Europe_west", "Europe_suture_zone", "Europe_east"]],
                    "k8_Americas_p1": [round(45.6885),
                                    ["Americas"]],
                    "k8_Americas_p2": [round(45.6885),
                                    ["mainland", "Caribbean"]],
                    "k8_mainland_p1": [round(49.7879),
                                       ["Americas"]],
                    "k8_mainland_p2": [round(49.7879),
                                       ["Americas", "United_States_east"]],
                    "k4_Transatlantic": [round(35.3440),
                                         ["Americas",
                                          "Europe_east", "Europe_west", 
                                          "Guinea", "Zambia"]],
                    "k4_Transatlantic_expandedAfr": [round(34.6543),
                                         ["Americas",
                                          "Europe_east", "Europe_west", 
                                          "Africa_west", "Africa_east"]],
                    "k4_Australia": [round(31.4914),
                                     ["Australia",
                                      "Europe_west", "Europe_east",
                                      "Africa_west", "Africa_east"]]}

def main():
    # Handle command line arguments
    popinfo_name_no_ext = sys.argv[1]

    # Set up config variables
    data_dir = "data/"
    popinfo_dir = f"{data_dir}popinfos/"
    sfs_dir = f"{data_dir}sfss/"
    vcf_file = f"{data_dir}vcfs/dest2_clustered_mainChroms.vcf.gz"
    n = 20

    construct_sfs(popinfo_name_no_ext, popinfo_dir, sfs_dir, vcf_file, n)

def construct_sfs(popinfo_name_no_ext, popinfo_dir, sfs_dir, vcf_file, n) -> None:
    """
    :param popinfo_name_no_ext str: Name of the popinfo file that lists samples
        to be incorporated into the constructed SFS and the populations to which
        they belong, without the file extension
    :param popinfo_dir str: Directory containing the popinfo file
    :param sfs_dir str: Directory in which to save the constructed SFS
    :param vcf_file str: Name of input VCF file, including path
    :param n int: Number of samples to which the SFS is projected down, i.e. the
        constructed SFS will have length `n+1` in each dimension. For this project,
        always chosen to be 20. This is distinct from effective sample size, which
        is calculated for each region (see `popinfo_metadata`) and used for rounding
        of non-integer allele counts that result from pool-seq.
    :return: None, but calls the function thatsaves the constructed SFS as a .npy 
        file in `sfs_dir`
    :rtype: None
    """
    #make_dirs([popinfo_dir, sfs_dir])
    mean_n_eff, pop_ids = popinfo_metadata[popinfo_name_no_ext]
    Npop = len(pop_ids)
    ns = [n] * Npop
    
    write_sfs(vcf_file=vcf_file, 
              popinfo_file=f"{popinfo_dir}{popinfo_name_no_ext}.popinfo",
              pop_ids=pop_ids,
              projections=ns,
              n_eff=mean_n_eff,
              rounding_method="counts_no_round_up",
              sfs_file=f"{sfs_dir}{popinfo_name_no_ext}.npy")

def write_sfs(vcf_file: str, 
              popinfo_file: str, 
              pop_ids: list, 
              projections: list, 
              n_eff: int, 
              rounding_method: str, 
              sfs_file: str) -> None:
    """
    :param vcf_file str: Name of input VCF file, including path
    :param popinfo_file str: Name of input popinfo file, including path
    :param pop_ids list: Ordered list of population names
    :param projections list: List of sample sizes at which to represent each population
        in the SFS. The SFS will be constructed such that its length equals `projections[i]+1`
        along the `i`th axis.
    :param n_eff int: Effective sample sizes, used in rounding of non-integer allele
        counts that result from pool-seq.
    :return: None, but saves the constructed SFS as a .npy file in `sfs_dir`
    :rtype: None
    """
    data_dict = make_data_dict_vcf_pooled(vcf_filename=vcf_file, 
                                          popinfo_filename=popinfo_file, 
                                          n_eff=n_eff,
                                          rounding_method=rounding_method)

    print("Data dictionary loaded. Now constructing site frequency spectrum...")
    fs = from_data_dict(data_dict, 
                        pop_ids=pop_ids,
                        projections=projections,
                        polarized=False)
    
    unused_pop_lists = [[2, 4], [1, 4], [2, 3], [1, 3]]
    for unused_pops in unused_pop_lists:
        print(fs.marginalize(unused_pops).pop_ids)
    
    print("SFS constructed. Now saving...")
    np.save(sfs_file, fs.data)

def from_data_dict(data_dict: dict, 
                   pop_ids: list, 
                   projections: list, 
                   mask_corners: bool=True,
                   polarized: bool=True) -> moments.Spectrum:
    """
    :param data_dict dict: A dictionary of SNP data containing all necessary information
        to construct a SFS. Format specified in moments documentation here: 
        https://moments.readthedocs.io/en/main/api/api_moments.html#moments.Spectrum.from_data_dict
    :param pop_ids list: Ordered list of population names
    :param projections list: List of sample sizes at which to represent each population
        in the SFS. The SFS will be constructed such that its shape equals `projections`,
        but with each element incremented to account for absent (i.e. multiplicity 
        0) alleles.
    :param mask_corners bool: Whether to mask the corners of the SFS
    :param polarized bool: If `False`, then fold the SFS because ancestral vs. derived
        allele states are unknown.
    :return: Folded moments.Spectrum object constructed from the `data_dict` and
        `pop_ids`, with size specified by `projections`
    :rtype: moments.Spectrum
    
    This function is a replicate of moments.from_data_dict, but stripped down to
    only what is needed for this project so that it can run significantly faster,
    especially in the absence of the overhead of importing moments.
    """

    Npops = len(pop_ids)
    fs = np.zeros(np.asarray(projections) + 1)
    l = len(data_dict.items())
    for snp, snp_info in data_dict.items():
        # Skip SNPs that aren't biallelic.
        if len(snp_info['segregating']) != 2:
            continue

        allele1, allele2 = snp_info['segregating']
        if not polarized:
            # If we don't want to polarize, we can choose which allele is
            # derived arbitrarily since we'll fold anyways.
            outgroup_allele = allele1
        elif 'outgroup_allele' in snp_info\
            and snp_info['outgroup_allele'] != '-'\
            and snp_info['outgroup_allele'] in snp_info['segregating']:
            # Otherwise we need to check that it's a useful outgroup
            outgroup_allele = snp_info['outgroup_allele']
        else: 
            # If we're polarized and we didn't have good outgroup info, skip
            # this SNP.
            continue

        # Extract the allele calls for each population.
        allele1_calls = [snp_info['calls'][pop][0] for pop in pop_ids]
        allele2_calls = [snp_info['calls'][pop][1] for pop in pop_ids]
        # How many chromosomes did we call successfully in each population?
        successful_calls = [a1 + a2 for (a1, a2) in zip(allele1_calls, allele2_calls)]

        # Which allele is derived (different from outgroup)?
        if allele1 == outgroup_allele:
            derived_calls = allele2_calls
        elif allele2 == outgroup_allele:
            derived_calls = allele1_calls

        # To handle arbitrary numbers of populations in the fs, we need
        # to do some tricky slicing.
        slices = [[np.newaxis] * len(pop_ids) for ii in range(Npops)]
        for ii in range(len(pop_ids)):
            slices[ii][ii] = slice(None,None,None)
    
        # Do the projection for this SNP.
        pop_contribs = []
        iter = zip(projections, successful_calls, derived_calls)
        for pop_ii, (p_to, p_from, hits) in enumerate(iter):
            contrib = moments.Numerics._cached_projection(p_to,p_from,hits)[tuple(slices[pop_ii])]
            pop_contribs.append(contrib)
        fs += functools.reduce(operator.mul, pop_contribs)
    fsout = moments.Spectrum(fs, mask_corners=mask_corners, 
                        pop_ids=pop_ids)
    if polarized:
        return fsout
    else:
        return fsout.fold()

# Adapted for pool-seq VCF files from moments.Misc.make_data_dict_vcf
def make_data_dict_vcf_pooled(vcf_filename: str, 
                              popinfo_filename: str,
                              n_eff: int, 
                              rounding_method: str='counts_no_round_up',
                              filter: bool=True, 
                              flanking_info: list=[None, None]) -> dict:
    """
    :param vcf_filename str: Name of input VCF file, including path
    :param popinfo_filename str: Name of input popinfo file, including path
    :param n_eff int: Effective sample sizes, used in rounding of non-integer allele 
        counts that result from pool-seq.
    :param rounding_method str: Method used to round non-integer allele counts.
        `counts_no_round_up` rounds to the nearest integer. `counts` rounds to the
        nearest integer, but always rounds up to 1 for allele counts in the interval
        (0, 1) in order to prevent "rounding out" rare alleles. `binom` performs
        binomial resampling, thus adding another layer of noise to the pool-seq 
        data.
    :param filter bool: Whether to filter out SNPs for which there is missing information
        in the VCF file
    :param flanking_info list: "Flanking information" to include about downstream 
        and upstream bases adjacent to SNPs
    :return: `data_dict` object that can be used as input to `moments.Spectrum.from_data_dict`,
        which I replicate here as `from_data_dict`.
    :rtype: dict

    This function adapts moments.Misc.make_data_dict_vcf to VCFs representing pool-seq 
    data. It is thus redundant with the function `dadi_inputs()` in the R package
    genomalicious (Thia and Riginos 2020, bioRxiv).
    """

    # Read population information from file based on extension
    if os.path.splitext(popinfo_filename)[1] == '.gz':
        import gzip
        popinfo_file = gzip.open(popinfo_filename)
    elif os.path.splitext(popinfo_filename)[1] == '.zip':
        import zipfile
        archive = zipfile.ZipFile(popinfo_filename)
        namelist = archive.namelist()
        if len(namelist) != 1:
            raise ValueError("Must be only a single popinfo file in zip "
                             "archive: {}".format(popinfo_filename))
        popinfo_file = archive.open(namelist[0])
    else:
        popinfo_file = open(popinfo_filename)
    # pop_dict has key, value pairs of "SAMPLE_NAME" : "POP_NAME"
    popinfo_dict = moments.Misc._get_popinfo(popinfo_file)
    popinfo_file.close()

    # Open VCF file
    if os.path.splitext(vcf_filename)[1] == '.gz':
        import gzip
        vcf_file = gzip.open(vcf_filename)
    elif os.path.splitext(vcf_filename)[1] == '.zip':
        import zipfile
        archive = zipfile.ZipFile(vcf_filename)
        namelist = archive.namelist()
        if len(namelist) != 1:
            raise ValueError("Must be only a single vcf file in zip "
                             "archive: {}".format(vcf_filename))
        vcf_file = archive.open(namelist[0])
    else:
        vcf_file = open(vcf_filename)
    
    data_dict = {}
    for line in vcf_file:
        # decoding lines for Python 3 - probably a better way to handle this
        try:
            line = line.decode()
        except AttributeError:
            pass
        # Skip metainformation
        if line.startswith('##'):
            continue
        # Read header
        if line.startswith('#'):
            header_cols = line.split()
            # Ensure there is at least one sample
            if len(header_cols) <= 9:
                raise ValueError("No samples in VCF file")
            # Use popinfo_dict to get the order of populations present in VCF
            poplist = [popinfo_dict[sample] if sample in popinfo_dict else None
                       for sample in header_cols[9:]] 
            continue
            
        # Read SNP data
        cols = line.split()
        snp_id = '_'.join(cols[:2]) # CHROM_POS
        snp_dict = {}
        
        # Skip SNP if filter is set to True and it fails a filter test
        if filter and cols[6] != 'PASS' and cols[6] != '.':
            continue

        # Add reference and alternate allele info to dict
        ref, alt = (allele.upper() for allele in cols[3:5])
        if ref not in ['A', 'C', 'G', 'T'] or alt not in ['A', 'C', 'G', 'T']:
            # Skip line if site is not an SNP
            continue
        snp_dict['segregating'] = (ref, alt)
        snp_dict['context'] = '-' + ref + '-'

        # Add ancestral allele information if available
        info = cols[7].split(';')
        for field in info:
            if field.startswith('AA'):
                outgroup_allele = field[3:].upper()
                if outgroup_allele not in ['A','C','G','T']:    
                    # Skip if ancestral not single base A, C, G, or T
                    outgroup_allele = '-'
                break
        else:
            outgroup_allele = '-'
        snp_dict['outgroup_allele'] = outgroup_allele
        snp_dict['outgroup_context'] = '-' + outgroup_allele + '-'

        # Add flanking info if it is present
        rflank, aflank = flanking_info
        for field in info:
            if rflank and field.startswith(rflank):
                flank = field[len(rflank+1):].upper()
                if not (len(flank) == 2 or len(flank) == 3):
                    continue
                prevb, nextb = flank[0], flank[-1]
                if prevb not in ['A','C','T','G']:
                    prevb = '-'
                if nextb not in ['A','C','T','G']:
                    nextb = '-'
                snp_dict['context'] = prevb + ref + nextb
                continue
            if aflank and field.startswith(aflank):
                flank = field[len(aflank+1):].upper()
                if not (len(flank) == 2 or len(flank) == 3):
                    continue
                prevb, nextb = flank[0], flank[-1]
                if prevb not in ['A','C','T','G']:
                    prevb = '-'
                if nextb not in ['A','C','T','G']:
                    nextb = '-'
                snp_dict['outgroup_context'] = prevb + outgroup_allele + nextb
        
        # Add reference and alternate allele calls for each population
        calls_dict = {}
        full_info = True
        gtindex = cols[8].split(':').index('GT')
        refindex = cols[8].split(':').index('RD')
        altindex = cols[8].split(':').index('AD')
        for pop, sample in zip(poplist, cols[9:]):
            if pop is None:
                continue
            gt = sample.split(':')[gtindex]
            g1, g2 = gt[0], gt[2]
            if g1 == '.' or g2 == '.':
                full_info = False
                break

            ref_count = int(sample.split(':')[refindex])
            alt_count = int(sample.split(':')[altindex])
            vaf = alt_count / (ref_count + alt_count)

            if rounding_method == 'counts_no_round_up':
                alt_count = np.round(vaf * n_eff)
            elif rounding_method == 'counts':
                unrounded_ac = vaf * n_eff
                alt_count = np.round(unrounded_ac) if unrounded_ac > 0.5 \
                       else 1 if unrounded_ac > 0 \
                       else 0
            elif rounding_method == 'binom':
                alt_count = np.random.binomial(n_eff, vaf)
            else:
                raise "Unrecognized rounding method."

            if pop not in calls_dict:
                calls_dict[pop] = (0,0)
            refcalls, altcalls = calls_dict[pop]
            refcalls += n_eff - alt_count
            altcalls += alt_count
            calls_dict[pop] = (refcalls, altcalls)

        if full_info:
            snp_dict['calls'] = calls_dict
            data_dict[snp_id] = snp_dict

    vcf_file.close()
    return data_dict

def make_dirs(dir_list: list) -> None:
    """
    :param dir_list list: List of directory paths to be created, if they do not 
        already exist
    
    :rtype: None

    Given a list of paths, make sure that those paths exist by creating non-existent
    terminal directories.
    """
    for dir in dir_list:
        if not os.path.exists(dir):
            os.makedirs(dir)

if __name__ == "__main__":
    main()
    
