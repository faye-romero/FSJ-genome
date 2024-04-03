#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="circos.GalGal"
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=50G
#SBATCH --time=3-00:00:00
#SBATCH --output=circos.GalGal.log

## Whole-genome circos plots ##
## Faye Romero. University of Rochester. 25 March 2024 ##

module load mummer/4.0.0beta1
module load gnuplot
module load samtools
module load minimap2

REF='GCA_024206055.2_GGswu_genomic_Chronly.fasta'
QUERY='AphCoe_V3_internal_Mar2024_Chronly.fasta'
PREFIX='GalGal_AphCoe'

# Make sure both REF and QUERY have .fai files. Get these from samtools faidx.

# minimap2
minimap2 -x asm20 $REF $QUERY > ${PREFIX}.aln.paf

# Create karyogram
#Ref has 40 scaffolds. 360/40 = 9. Therefore, set h+= to 9.
awk '{hue=sprintf("%03d", h); print "chr - "$1" "$1" 0 "$2" hue"hue; h+=9;}' ${REF}.fai >  ${PREFIX}.chr.kar
awk '{print "chr - "$1" "$1" 0 "$2" white"}' ${QUERY}.fai >> ${PREFIX}.chr.kar

# Create links from a minimap2 aln
awk 'BEGIN {FS=OFS="\t"} {print $0, $4 - $3}' ${PREFIX}.aln.paf | sort -k19,19nr | awk 'BEGIN {OFS="\t"; z=0} {z+=1; print $6, $8, $9, $1, $3, $4, "z="z}' > ${PREFIX}.links.tsv

# Order the karyogram file
awk '{hue=sprintf("%03d", h); print "chr - "$1" "$1" 0 "$2" hue"hue; h+=9;}' ${REF}.fai >  ${PREFIX}.chr.kar
tac ${QUERY}.fai | awk '{print "chr - "$1" "$1" 0 "$2" white"}' >> ${PREFIX}.chr.kar

# Add color
perl addCol.pl ${PREFIX}.chr.kar ${PREFIX}.links.tsv > ${PREFIX}.newlinks.tsv

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
label_size       = 30p
</ideogram>

# LINKS
<links>
<link>
file                = ${PREFIX}.newlinks.tsv
thickness           = 40p
radius              = 0.95r
ribbon              = yes
flat		    	= yes
bezier_radius        = 0r
bezier_radius_purity = 0.5
color               = white
stroke_color        = black
stroke_thickness    = 0.02p
</link>
</links>"

echo "$CONFIG" > ${PREFIX}.circos.conf

# Then, save chr.kar, newlinks.tsv, and circos.conf to your local computer.