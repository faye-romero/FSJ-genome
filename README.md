# Florida Scrub-Jay genome assembly and annotation
  
Scripts for *de novo* genome assembly of the Florida Scrub-Jay (*Aphelocoma coerulescens*) genome.  
  
Find the genome assembly under NCBI accession **GCA_041296385.1** (UR_Acoe_1.0). Raw sequence reads are available on NCBI under BioProjects PRJNA1076903 and PRJNA1097984. The genome annotation, repeat library and repeat annotation, and sex-averaged and sex-specific linkage maps are available on Figshare at figshare(.)com/projects/Florida_Scrub-Jay_genome_assembly/220939.  
  
If you use any data associated with this project, please cite:  
  
Faye G. Romero, Felix E.G. Beaudry, Eyvind Hovmand Warner, Tram N. Nguyen, John W. Fitzpatrick, Nancy Chen. **A new high-quality genome assembly and annotation for the threatened Florida Scrub-Jay** (*Aphelocoma coerulescens*). *bioRxiv* 2024.04.05.588142. DOI: 10.1101/2024.04.05.588142.  
  
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
  
Data in this repository (found in `data/`) include:  
- `TE_lengths_v123_July2024.Rdata`: repeat element classifications and lengths for 3 versions of the Florida Scrub-Jay genome (v1 = Feng et al. 2020, Nature; v2 = Driscoll & Beaudry et al. 2024, Genetics; v3 = Romero et al. 2024, bioRxiv)  
- `TE_counts_v23_July2024.Rdata`: repeat element classifications and counts for 2 versions of the Florida Scrub-Jay genome (v2 = Driscoll & Beaudry et al. 2024, Genetics; v3 = Romero et al. 2024, bioRxiv)  
- `FSJgenome_July2024_FINAl_*.fasta.tbl`: RepeatMasker table summaries  
- `AllSamples_AvgCov_Wunloc3_AllPos_1kbwin.Rdata`: average read depth for 50 jays, in 1 kb windows, across the chimeric ZW contig  
- `Indiv_Cov_Per_Scaf.Rdata`: average read depth for 50 jays across all scaffolds in the genome  
- `unmasked_seq_pos_FSJgenome_July2024_FINAL.bed`: positions of unmasked (lowercase) sequence in the genome (generated with the script `sex-linked-scaffold-id/extract_unmasked_seq_pos.py`)  
  
---
  
Please direct any questions about this project to Faye Romero (fromero3@ur.rochester.edu).
