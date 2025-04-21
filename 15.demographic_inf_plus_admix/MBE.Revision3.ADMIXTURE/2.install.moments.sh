module load python3.12-anaconda/2024.06-1


conda create \
-n moments_jcbn \
python=3.12 \
dadi \
ipykernel \
-c conda-forge

# Activates conda environemnt
conda activate moments_jcbn

#Installs moments
conda install moments -c bioconda

# install pandas
conda install pandas