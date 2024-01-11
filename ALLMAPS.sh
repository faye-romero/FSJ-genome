## Steps to run ALLMAPS linkage map-based scaffolder on a genome assembly ##
## Faye Romero, University of Rochester ##
## 25 Sept 2023 ##

###### Ran the following on command line via an interactive session ######
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=200G

# Install instructions: github(dot)com/tanghaibao/jcvi/wiki/ALLMAPS%3A-How-to-install
# I chose to install Python dependencies in a conda environment

# Load modules + conda environment and set PATH
module unload miniconda miniconda3 python python3
module load miniconda3
source activate allmaps # Activate conda environment with Python dependencies
export PATH=/home/fromero3/.local/lib/python3.9/site-packages/jcvi:$PATH # ALLMAPS path
export PATH=/software/liftover/220922:$PATH # liftOver path
module load texlive

# Set wd
cd /scratch/nchen11_lab/Faye/ReferenceAug23/frankenstein_H1172

# Set variables
MAP='/scratch/nchen11_lab/Faye/ReferenceAug23/beadchip_out/Hifiasm_dip.hap1.H1172pre_salsafrank_beadchipSeqLoc_allmap.csv' # linkage map
GENOME='/scratch/nchen11_lab/Faye/ReferenceAug23/frankenstein_H1172/Hifiasm_dip.hap1.H1172pre_salsafrank.fasta' # genome
PREFIX='Hifiasm_dip.hap1.H1172pre_salsafrank' # name of assembly

# Prepare ALLMAPS input
python3 -m jcvi.assembly.allmaps merge ${MAP} -o ${PREFIX}_MAP.bed -w weights.txt
# Note that you need to name the .bed file something like "genome_MAP.bed", and not just the genome name, to avoid weird errors

# Run ALLMAPS
python3 -m jcvi.assembly.allmaps path ${PREFIX}_MAP.bed ${GENOME} -w weights.txt --cpus=4

source deactivate allmaps
