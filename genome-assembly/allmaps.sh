#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="FSJv3.allmaps"
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=200G
#SBATCH --time=3-00:00:00
#SBATCH --output=allmaps.03042024.log

## ALLMAPS ##
## Faye Romero. University of Rochester. 04 March 2024 ##

# Run on command line via an interactive session

# Load modules and conda envs
module unload miniconda miniconda3 python python3
module load miniconda3
source activate allmaps # Activate conda env
export PATH=/home/fromero3/.local/lib/python3.9/site-packages/jcvi:$PATH #Path to JCVI utils
export PATH=/software/liftover/220922:$PATH #liftOver path
module load texlive

echo "$(date).....begin ALLMAPS"

# Set wd and variables
cd /scratch/nchen11_lab/new_FSJgenome/allmaps_out_FIN
MAP='/scratch/nchen11_lab/new_FSJgenome/allmaps_out_FIN/AphCoe_V3_beadchipSeqLoc_allmap.csv'
GENOME='/scratch/nchen11_lab/new_FSJgenome/allmaps_out_FIN/AphCoe_V3.fasta'
PREFIX='AphCoe_V3'

# Prepare ALLMAPS input
python3 -m jcvi.assembly.allmaps merge ${MAP} -o ${PREFIX}_MAP.bed -w weights.txt
#Note that you need to name the .bed file something like "genome_MAP.bed", and not just the genome name, to avoid weird errors!

#Run ALLMAPS
python3 -m jcvi.assembly.allmaps path ${PREFIX}_MAP.bed ${GENOME} -w weights.txt --cpus=4

source deactivate
echo "$(date).....end ALLMAPS"

######################################################################

# 'allmaps' conda env specs below
# packages in environment at /home/fromero3/.conda/envs/allmaps:
#
# Name                    Version                   Build  Channel
# _libgcc_mutex             0.1                 conda_forge    conda-forge
# _openmp_mutex             4.5                       2_gnu    conda-forge
# absl-py                   1.4.0                    pypi_0    pypi
# alsa-lib                  1.2.8                h166bdaf_0    conda-forge
# attr                      2.5.1                h166bdaf_1    conda-forge
# biopython                 1.81            py310h1fa729e_0    conda-forge
# brewer2mpl                1.4.1                      py_3    conda-forge
# brotli                    1.0.9                h166bdaf_9    conda-forge
# brotli-bin                1.0.9                h166bdaf_9    conda-forge
# brotli-python             1.0.9           py310hd8f1fbe_9    conda-forge
# bx-python                 0.10.0                   pypi_0    pypi
# bzip2                     1.0.8                h7f98852_4    conda-forge
# ca-certificates           2023.7.22            hbcca054_0    conda-forge
# cairo                     1.16.0            hbbf8b49_1016    conda-forge
# certifi                   2023.7.22          pyhd8ed1ab_0    conda-forge
# charset-normalizer        3.2.0              pyhd8ed1ab_0    conda-forge
# contourpy                 1.1.0           py310hd41b1e2_0    conda-forge
# crossmap                  0.6.6                    pypi_0    pypi
# cycler                    0.11.0             pyhd8ed1ab_0    conda-forge
# cython                    3.0.2                    pypi_0    pypi
# dbus                      1.13.6               h5008d03_3    conda-forge
# deap                      1.4.1           py310h7cbd5c2_0    conda-forge
# expat                     2.5.0                hcb278e6_1    conda-forge
# font-ttf-dejavu-sans-mono 2.37                 hab24e00_0    conda-forge
# font-ttf-inconsolata      3.000                h77eed37_0    conda-forge
# font-ttf-source-code-pro  2.038                h77eed37_0    conda-forge
# font-ttf-ubuntu           0.83                 hab24e00_0    conda-forge
# fontconfig                2.14.2               h14ed4e7_0    conda-forge
# fonts-conda-ecosystem     1                             0    conda-forge
# fonts-conda-forge         1                             0    conda-forge
# fonttools                 4.42.1          py310h2372a71_0    conda-forge
# freetype                  2.12.1               hca18f0e_1    conda-forge
# ftpretty                  0.4.0                    pypi_0    pypi
# gettext                   0.21.1               h27087fc_0    conda-forge
# glib                      2.76.4               hfc55251_0    conda-forge
# glib-tools                2.76.4               hfc55251_0    conda-forge
# graphite2                 1.3.13            h58526e2_1001    conda-forge
# gst-plugins-base          1.22.3               h938bd60_1    conda-forge
# gstreamer                 1.22.3               h977cf35_1    conda-forge
# harfbuzz                  7.3.0                hdb3a94d_0    conda-forge
# icu                       72.1                 hcb278e6_0    conda-forge
# idna                      3.4                pyhd8ed1ab_0    conda-forge
# imageio                   2.31.2                   pypi_0    pypi
# jcvi                      1.3.7                    pypi_0    pypi
# keyutils                  1.6.1                h166bdaf_0    conda-forge
# kiwisolver                1.4.5           py310hd41b1e2_0    conda-forge
# krb5                      1.20.1               h81ceb04_0    conda-forge
# lame                      3.100             h166bdaf_1003    conda-forge
# lazy-loader               0.3                      pypi_0    pypi
# lcms2                     2.15                 haa2dc70_1    conda-forge
# ld_impl_linux-64          2.40                 h41732ed_0    conda-forge
# lerc                      4.0.0                h27087fc_0    conda-forge
# libblas                   3.9.0           17_linux64_openblas    conda-forge
# libbrotlicommon           1.0.9                h166bdaf_9    conda-forge
# libbrotlidec              1.0.9                h166bdaf_9    conda-forge
# libbrotlienc              1.0.9                h166bdaf_9    conda-forge
# libcap                    2.69                 h0f662aa_0    conda-forge
# libcblas                  3.9.0           17_linux64_openblas    conda-forge
# libclang                  16.0.6          default_h1cdf331_1    conda-forge
# libclang13                16.0.6          default_h4d60ac6_1    conda-forge
# libcups                   2.3.3                h36d4200_3    conda-forge
# libdeflate                1.18                 h0b41bf4_0    conda-forge
# libedit                   3.1.20191231         he28a2e2_2    conda-forge
# libevent                  2.1.12               hf998b51_1    conda-forge
# libexpat                  2.5.0                hcb278e6_1    conda-forge
# libffi                    3.4.2                h7f98852_5    conda-forge
# libflac                   1.4.3                h59595ed_0    conda-forge
# libgcc-ng                 13.1.0               he5830b7_0    conda-forge
# libgcrypt                 1.10.1               h166bdaf_0    conda-forge
# libgfortran-ng            13.1.0               h69a702a_0    conda-forge
# libgfortran5              13.1.0               h15d22d2_0    conda-forge
# libglib                   2.76.4               hebfc3b9_0    conda-forge
# libgomp                   13.1.0               he5830b7_0    conda-forge
# libgpg-error              1.47                 h71f35ed_0    conda-forge
# libiconv                  1.17                 h166bdaf_0    conda-forge
# libjpeg-turbo             2.1.5.1              h0b41bf4_0    conda-forge
# liblapack                 3.9.0           17_linux64_openblas    conda-forge
# libllvm16                 16.0.6               h5cf9203_2    conda-forge
# libnsl                    2.0.0                h7f98852_0    conda-forge
# libogg                    1.3.4                h7f98852_1    conda-forge
# libopenblas               0.3.23          pthreads_h80387f5_0    conda-forge
# libopus                   1.3.1                h7f98852_1    conda-forge
# libpng                    1.6.39               h753d276_0    conda-forge
# libpq                     15.3                 hbcd7760_1    conda-forge
# libsndfile                1.2.2                hbc2eb40_0    conda-forge
# libsqlite                 3.43.0               h2797004_0    conda-forge
# libstdcxx-ng              13.1.0               hfd8a6a1_0    conda-forge
# libsystemd0               254                  h3516f8a_0    conda-forge
# libtiff                   4.5.1                h8b53f26_1    conda-forge
# libuuid                   2.38.1               h0b41bf4_0    conda-forge
# libvorbis                 1.3.7                h9c3ff4c_0    conda-forge
# libwebp-base              1.3.1                hd590300_0    conda-forge
# libxcb                    1.15                 h0b41bf4_0    conda-forge
# libxkbcommon              1.5.0                h5d7e998_3    conda-forge
# libxml2                   2.11.5               h0d562d8_0    conda-forge
# libzlib                   1.2.13               hd590300_5    conda-forge
# lz4-c                     1.9.4                hcb278e6_0    conda-forge
# markdown-it-py            3.0.0              pyhd8ed1ab_0    conda-forge
# matplotlib                3.7.2           py310hff52083_0    conda-forge
# matplotlib-base           3.7.2           py310hf38f957_0    conda-forge
# mdurl                     0.1.0              pyhd8ed1ab_0    conda-forge
# more-itertools            10.1.0             pyhd8ed1ab_0    conda-forge
# mpg123                    1.31.3               hcb278e6_0    conda-forge
# munkres                   1.1.4              pyh9f0ad1d_0    conda-forge
# mysql-common              8.0.33               hf1915f5_2    conda-forge
# mysql-libs                8.0.33               hca2cd23_2    conda-forge
# natsort                   8.4.0              pyhd8ed1ab_0    conda-forge
# ncurses                   6.4                  hcb278e6_0    conda-forge
# networkx                  3.1                pyhd8ed1ab_0    conda-forge
# nspr                      4.35                 h27087fc_0    conda-forge
# nss                       3.92                 h1d7d5a4_0    conda-forge
# numpy                     1.25.2          py310ha4c1d20_0    conda-forge
# openjpeg                  2.5.0                hfec8fc6_2    conda-forge
# openssl                   3.1.2                hd590300_0    conda-forge
# ortools                   9.7.2996                 pypi_0    pypi
# packaging                 23.1               pyhd8ed1ab_0    conda-forge
# pcre2                     10.40                hc3806b6_0    conda-forge
# pillow                    10.0.0          py310h582fbeb_0    conda-forge
# pip                       23.2.1             pyhd8ed1ab_0    conda-forge
# pixman                    0.40.0               h36c2ea0_0    conda-forge
# platformdirs              3.10.0             pyhd8ed1ab_0    conda-forge
# ply                       3.11                       py_1    conda-forge
# pooch                     1.7.0              pyha770c72_3    conda-forge
# protobuf                  4.24.2                   pypi_0    pypi
# pthread-stubs             0.4               h36c2ea0_1001    conda-forge
# pulseaudio-client         16.1                 hb77b528_4    conda-forge
# pybedtools                0.9.1                    pypi_0    pypi
# pybigwig                  0.3.22                   pypi_0    pypi
# pygments                  2.16.1             pyhd8ed1ab_0    conda-forge
# pyparsing                 3.0.9              pyhd8ed1ab_0    conda-forge
# pyqt                      5.15.9          py310h04931ad_4    conda-forge
# pyqt5-sip                 12.12.2         py310hc6cd4ac_4    conda-forge
# pysam                     0.21.0                   pypi_0    pypi
# pysocks                   1.7.1              pyha2e5f31_6    conda-forge
# python                    3.10.12         hd12c33a_0_cpython    conda-forge
# python-dateutil           2.8.2              pyhd8ed1ab_0    conda-forge
# python-graphviz           0.20.1                   pypi_0    pypi
# python_abi                3.10                    3_cp310    conda-forge
# pywavelets                1.4.1                    pypi_0    pypi
# qt-main                   5.15.8              h01ceb2d_12    conda-forge
# readline                  8.2                  h8228510_1    conda-forge
# requests                  2.31.0             pyhd8ed1ab_0    conda-forge
# rich                      13.5.1             pyhd8ed1ab_0    conda-forge
# scikit-image              0.21.0                   pypi_0    pypi
# scipy                     1.11.2          py310ha4c1d20_0    conda-forge
# setuptools                68.1.2             pyhd8ed1ab_0    conda-forge
# sip                       6.7.11          py310hc6cd4ac_0    conda-forge
# six                       1.16.0             pyh6c4a22f_0    conda-forge
# tifffile                  2023.8.25                pypi_0    pypi
# tk                        8.6.12               h27826a3_0    conda-forge
# toml                      0.10.2             pyhd8ed1ab_0    conda-forge
# tomli                     2.0.1              pyhd8ed1ab_0    conda-forge
# tornado                   6.3.3           py310h2372a71_0    conda-forge
# typing-extensions         4.7.1                hd8ed1ab_0    conda-forge
# typing_extensions         4.7.1              pyha770c72_0    conda-forge
# tzdata                    2023c                h71feb2d_0    conda-forge
# unicodedata2              15.0.0          py310h5764c6d_0    conda-forge
# urllib3                   2.0.4              pyhd8ed1ab_0    conda-forge
# webcolors                 1.13                     pypi_0    pypi
# wheel                     0.41.2             pyhd8ed1ab_0    conda-forge
# xcb-util                  0.4.0                hd590300_1    conda-forge
# xcb-util-image            0.4.0                h8ee46fc_1    conda-forge
# xcb-util-keysyms          0.4.0                h8ee46fc_1    conda-forge
# xcb-util-renderutil       0.3.9                hd590300_1    conda-forge
# xcb-util-wm               0.4.1                h8ee46fc_1    conda-forge
# xkeyboard-config          2.39                 hd590300_0    conda-forge
# xorg-kbproto              1.0.7             h7f98852_1002    conda-forge
# xorg-libice               1.1.1                hd590300_0    conda-forge
# xorg-libsm                1.2.4                h7391055_0    conda-forge
# xorg-libx11               1.8.6                h8ee46fc_0    conda-forge
# xorg-libxau               1.0.11               hd590300_0    conda-forge
# xorg-libxdmcp             1.1.3                h7f98852_0    conda-forge
# xorg-libxext              1.3.4                h0b41bf4_2    conda-forge
# xorg-libxrender           0.9.11               hd590300_0    conda-forge
# xorg-renderproto          0.11.1            h7f98852_1002    conda-forge
# xorg-xextproto            7.3.0             h0b41bf4_1003    conda-forge
# xorg-xf86vidmodeproto     2.3.1             h7f98852_1002    conda-forge
# xorg-xproto               7.0.31            h7f98852_1007    conda-forge
# xz                        5.2.6                h166bdaf_0    conda-forge
# zlib                      1.2.13               hd590300_5    conda-forge
# zstd                      1.5.5                hfc55251_0    conda-forge