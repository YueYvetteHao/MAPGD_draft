#!/bin/bash
#SBATCH --job-name=mapgd	    # Job name
#SBATCH --ntasks=1                  # CPU
#SBATCH --mem=5gb                   # Job memory request
#SBATCH --time=24:00:00             # Time limit hrs:min:sec
#SBATCH --output=mapgd_%j.log       # Standard output and error log


# Load necessary modules
module load GSL/2.6-GCC-9.3.0
module load samtools/1.9

# Working directory
cd /scratch/yhao/MAPGD_plants

# Download input bam files
wget https://1001genomes.org/data/JGI/JGIHeazlewood2011/releases/current/TAIR10/strains/Alc-0/Alc-0.bam --no-check-certificate
wget https://1001genomes.org/data/JGI/JGIHeazlewood2011/releases/current/TAIR10/strains/Jea/Jea.bam --no-check-certificate
wget https://1001genomes.org/data/JGI/JGIHeazlewood2011/releases/current/TAIR10/strains/Oy-0/Oy-0.bam --no-check-certificate
wget https://1001genomes.org/data/JGI/JGIHeazlewood2011/releases/current/TAIR10/strains/Ri-0/Ri-0.bam --no-check-certificate
wget https://1001genomes.org/data/JGI/JGIHeazlewood2011/releases/current/TAIR10/strains/Sakata/Sakata.bam --no-check-certificate

# Input file list
ls *.bam > filelist.txt
sed -i 's/.bam//' filelist.txt

# Sort bam files and filter for properly mapped reads
input='filelist.txt'
while IFS= read -r line
do
        echo $line
        samtools sort -o $line.sort.bam $line.bam
        # SAM flags: https://broadinstitute.github.io/picard/explain-flags.html
        samtools view -q 20 -f 3 -F 3844 -b $line.sort.bam > $line.filtered.sort.bam
        samtools index $line.filtered.sort.bam
done < $input

# Get header file
samtools view -H Alc-0.sort.bam > Arabidopsis.header

# Generate mpileup file
samtools mpileup -q 25 -Q 25 -B Alc-0.filtered.sort.bam Jea.filtered.sort.bam Oy-0.filtered.sort.bam Ri-0.filtered.sort.bam Sakata.filtered.sort.bam > Arabidopsis.mpileup

# Convert .mpileup file to .pro file
mapgd proview -i Arabidopsis.mpileup -H Arabidopsis.header > Arabidopsis.pro
