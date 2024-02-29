import os

"""
This script outputs an options file for a Slurm array job to run via `wrapper_construct_jackknife_sfs.slurm`.
In its present state, it is set up to only write the options file for jackknifed
replicates of the k4_Transatlantic region.
"""



def main():
    # Setup
    data_dir = "/project/berglandlab/david_DEST2/"
    popinfo_dir = data_dir + "popinfos/"
    options_file = f"options_construct_jackknife_sfs.txt"

    output_string = ""

    for jk_popinfo_dir in os.listdir(popinfo_dir):
        # Filter out non-directory files
        if "." in jk_popinfo_dir:
            continue

        # Prevent repetition of previous array jobs from this options file
        if jk_popinfo_dir not in ["k4_Transatlantic_jackknife"]:
            continue

        for popinfo_name in os.listdir(popinfo_dir + jk_popinfo_dir):
            if popinfo_name.endswith(".popinfo"):
                popinfo_name_no_ext = popinfo_name.rsplit(".", 1)[0]
                output_string += popinfo_name_no_ext + "\n"

    with open(options_file, 'w') as f:
        f.write(output_string.strip())

if __name__ == "__main__":
    main()
