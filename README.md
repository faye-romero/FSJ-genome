# Florida Scrub-Jay genome assembly and annotation
  
Scripts for *de novo* genome assembly of the Florida Scrub-Jay (*Aphelocoma coerulescens*) genome.  
  
Find the genome assembly under NCBI accession GCA_041296385.1, BioProject PRJNA1076903.  
  
If you use any data associated with this project, please cite:  
  
Faye G. Romero, Felix E.G. Beaudry, Eyvind Hovmand Warner, Tram N. Nguyen, John W. Fitzpatrick, Nancy Chen. **A new high-quality genome assembly and annotation for the threatened Florida Scrub-Jay** (*Aphelocoma coerulescens*). *bioRxiv* 2024.04.05.588142; doi: https(:)//doi(.)org/10.1101/2024.04.05.588142.  
  
---
  
Input data for basic genome assembly and annotation:  
* PacBio HiFi raw reads  
* Illumina paired-end reads of parents of the reference individual  
* Hi-C reads  
* linkage map  
* RNA-seq reads  
  
Input data for sex-linked scaffold identification pipeline:  
* Paired-end whole-genome Illumina reads for at least 1 male and 1 female individual (here, I use 25 male and 25 female)  
  
---
  
Other data in this repository include:  
* The genome annotation (`gene-annotation/UR_Acoe_1.0_genomic_annotation.gff3.gz`)  
* The master repeat library used for RepeatMasker annotation (`repeat-annotation/FSJgenome_July2024_FINAL_master_repeats.lib`)  
* Sex-averaged (`figures/FigureS1/FileS1.tsv`) and sex-specific (`figures/Figure2_LinkageMap/LOD3LOD5_anon_map.csv`) linkage maps for the Florida Scrub-Jay  
* Average coverage per scaffold for 50 Florida Scrub-Jays (`figures/Figure4_SexLinkedScaffs/Indiv_Cov_Per_Scaf.Rdata`)  
* Average coverage in 1 kb windows across a chimeric contig (`chimera-discovery/AllSamples_AvgCov_Wunloc3_AllPos_1kbwin.Rdata`)  
* RepeatMasker outputs for this genome assembly (`figures/Figure5_Repeats` and `figures/Supp_Figures/FigureS6`)  
  
Please direct any questions about this project to Faye Romero (fromero3@ur.rochester.edu).
