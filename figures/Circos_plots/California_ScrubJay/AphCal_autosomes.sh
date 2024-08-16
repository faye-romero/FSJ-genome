#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="circos.aphcal"
#SBATCH -N 1
#SBATCH -c 3
#SBATCH --nodelist=bhd0056
#SBATCH --mem=30G
#SBATCH --time=3-00:00:00
#SBATCH --output=circos.aphcal.log

# Interspecies whole-genome circos plots
# Adapted from https(:)//bioinf(.)cc/misc/2020/08/08/circos-ribbons(.)html

module load mummer/4.0.0beta1
module load gnuplot
module load samtools
module load minimap2

cd /scratch/nchen11_lab/Faye/new_FSJgenome

REF='/scratch/nchen11_lab/Faye/new_FSJgenome/autosomes_ORDERED/AphCal_autosomes_ORDERED.fasta'
QUERY='/scratch/nchen11_lab/Faye/new_FSJgenome/autosomes_ORDERED/FSJ_autosomes_ORDERED.fasta'
PREFIX='AphCal_AphCoe_autosomes'

mkdir -p ./circos_plots/${PREFIX}
cd ./circos_plots/${PREFIX}

# minimap2
minimap2 -x asm20 -t 3 $REF $QUERY > ${PREFIX}.aln.paf

# Create karyogram
#Ref has 31 scaffolds. 360/31 = 11.6. Therefore, set h+= to 12.
awk '{hue=sprintf("%03d", h); print "chr - "$1" "$1" 0 "$2" hue"hue; h+=12;}' ${REF}.fai >  ${PREFIX}.chr.kar
awk '{print "chr - "$1" "$1" 0 "$2" white"}' ${QUERY}.fai >> ${PREFIX}.chr.kar

# Create links from a minimap2 aln
awk 'BEGIN {FS=OFS="\t"} {print $0, $4 - $3}' ${PREFIX}.aln.paf | sort -k19,19nr | awk 'BEGIN {OFS="\t"; z=0} {z+=1; print $6, $8, $9, $1, $3, $4, "z="z}' > ${PREFIX}.links.tsv

# Order the karyogram file
awk '{hue=sprintf("%03d", h); print "chr - "$1" "$1" 0 "$2" hue"hue; h+=12;}' ${REF}.fai >  ${PREFIX}.chr.kar
tac ${QUERY}.fai | awk '{print "chr - "$1" "$1" 0 "$2" white"}' >> ${PREFIX}.chr.kar

# Add color
perl /scratch/nchen11_lab/Faye/new_FSJgenome/addCol.pl ${PREFIX}.chr.kar ${PREFIX}.links.tsv > ${PREFIX}.newlinks.tsv

# Finally, create config file circos.conf.

CONFIG="karyotype = ${PREFIX}.chr.kar
chromosomes_units = 100000
chromosomes_reverse = /Chr/

<<include colors_fonts_patterns.conf>>
<<include housekeeping.conf>>

# IMAGE
<image>
<<include image.conf>>
</image>

# IDEOGRAM
<ideogram>
<spacing>
default = 0.5u
break = 1u
</spacing>
radius           = 0.75r
thickness        = 100p
fill             = yes
stroke_color     = black
stroke_thickness = 2p
show_label       = yes
label_radius     = 1.05r
label_size       = 40p
</ideogram>

# LINKS
<links>
<link>
file                = ${PREFIX}.800kb.newlinks.tsv
thickness           = 40p
radius              = 0.95r
ribbon              = yes
flat		    	= yes
bezier_radius        = 0r
bezier_radius_purity = 0.5
color               = white
stroke_color        = black
stroke_thickness    = 1p
</link>
</links>"

echo "$CONFIG" > ${PREFIX}.circos.conf

# Then, save chr.kar, newlinks.tsv, and circos.conf to your local computer.