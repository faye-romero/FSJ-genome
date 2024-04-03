#!/bin/bash
#SBATCH -p rosalind
#SBATCH --job-name="hifiasm.triobin"
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=200G
#SBATCH --time=2-00:00:00
#SBATCH --output=hifiasm_triobin.08182023.log

## Trio-based hifiasm ##
## Faye Romero. University of Rochester. 18 August 2023 ##

#Note: To do trio-based genome assembly, you need yak. I installed it manually in /scratch/nchen11_lab/yak.

# Load modules
module load hifiasm/0.16.1
export PATH=$PATH:/scratch/nchen11_lab/yak #path to yak

# Set/make working directories
cd /scratch/nchen11_lab/Faye/ReferenceMarch23
#mkdir ./hifiasm_triobin_out_03202023

# Set variables
REFFQ=/scratch/nchen11_lab/FSJgenome/HiFi/FSJ.hifi.cut.fastq.gz #Ref raw reads
DADFQ=/scratch/nchen11_lab/Faye/ReferenceMarch23/pear_out/cleaned.merged.H1172.assembled.fastq #Dad raw reads
MOMFQ=/scratch/nchen11_lab/Faye/ReferenceMarch23/pear_out/cleaned.merged.H914.assembled.fastq #Mom raw reads

##count k-mers first with yak first:
/scratch/nchen11_lab/yak/yak count -k31 -b37 -t16 -o H1172.paternal.yak $DADFQ
/scratch/nchen11_lab/yak/yak count -k31 -b37 -t16 -o H914.maternal.yak $MOMFQ
echo "k-mer counting done!"
echo $(date)

##run primary assembly
echo $(date)
hifiasm -o ./hifiasm_triobin_out_03202023/FSJ_hifiasm_triobin_03202023.asm --primary -t 32 $REFFQ 2> ./hifiasm_triobin_out_03202023/FSJ_hifiasm_triobin_03202023.asm.pri.log
echo "primary hifiasm assembly done!"
echo $(date)

##run trio-based assembly
echo $(date)
hifiasm -o ./hifiasm_triobin_out_03202023/FSJ_hifiasm_triobin_03202023.asm -t 32 -1 ./H1172.paternal.yak -2 ./H914.maternal.yak /dev/null 2> ./hifiasm_triobin_out_03202023/FSJ_hifiasm_triobin_03202023.asm.trio.log
echo "trio binned hifiasm assemblies done!"
echo $(date)


## The phasing switch error rate and hamming error rate are able to be evaluated quickly by yak

cd /scratch/nchen11_lab/Faye/ReferenceMarch23/hifiasm_triobin_out_03202023

# Set variables
MAT='/scratch/nchen11_lab/Faye/ReferenceMarch23/H914.maternal.yak'
PAT='/scratch/nchen11_lab/Faye/ReferenceMarch23/H1172.paternal.yak'
PATHAP1='/scratch/nchen11_lab/Faye/ReferenceMarch23/hifiasm_triobin_out_03202023/FSJ_hifiasm_triobin_03202023.asm.dip.hap1.p_ctg.gfa'
MATHAP2='/scratch/nchen11_lab/Faye/ReferenceMarch23/hifiasm_triobin_out_03202023/FSJ_hifiasm_triobin_03202023.asm.dip.hap2.p_ctg.gfa'

# Combine paternal and maternal haplotype assemblies into a single fasta
cat $PATHAP1 |awk '{if (match($1, "^S")) { print ">"$2; print $3}}'|fold -c > trio.combined.fasta
cat $MATHAP2 |awk '{if (match($1, "^S")) { print ">"$2; print $3}}'|fold -c >> trio.combined.fasta

# Calculate switch/hamming error rate using yak
/scratch/nchen11_lab/yak/yak trioeval -t16 -c 5 -d 25 $PAT $MAT ./trio.combined.fasta | tail -n 3 > trioeval.out.yak.txt