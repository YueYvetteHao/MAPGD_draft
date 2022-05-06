#!/bin/bash
#SBATCH --job-name=mapgd	    # Job name
#SBATCH --ntasks=4                  # CPU
#SBATCH --mem=10gb                   # Job memory request
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

# Make a pro file of nucleotide-read quartets (counts of A, C, G, and T) from the mpileup file
mapgd proview -i Arabidopsis.mpileup -H Arabidopsis.header > Arabidopsis.pro

# Run the allele command to estimate allele and genotype frequencies from the pro file.
mapgd allele -i Arabidopsis.pro -o Arabidopsis -p Arabidopsis.clean
# .map file headline:
#@SCFNAME    	POS	REF	MAJOR	MINOR	COVERAG	MJ_FREQ	VR_FREQ	ERROR	NULL_ER	NULL_E2	F_STAT	MM_FREQ	Mm_FREQ	mm_FREQ	HETERO	POLY_LR	HWE_LR	GOF	EF_CHRM	IND_INC	IND_CUT	BEST_LL

# Determine the minimum and maximum population-coverage cut-off values by making a histogram of the population coverage.
awk '{print $6}' Arabidopsis.map > Arabidopsis_coverage.txt
# Make the coverage histogram in R.

# Run the filter command to filter the map file of ML estimates of the parameters.
mapgd filter -i Arabidopsis.map -p 20 -q 0.05 -Q 0.45 -c 10 -C 250 -o filtered_Arabidopsis
#-p: minimum value of the likelihood-ratio test statistic for polymorphism 
#-q: minimum minor-allele frequency estimate
#-Q: maximum minor-allele frequency estimate
#-c: minimum population coverage
#-C: maximum population coverage

# Run the genotype command to generate a file of genotype likelihoods.
mapgd genotype -p Arabidopsis.clean.pro -m filtered_Arabidopsis.map > Arabidopsis.genotype

# Format genotype file for relatedness analysis
# Remove the unnecessary header and footer
awk '{if ($3 != "MN_FREQ" && $3 >= 0.0 && $3 <= 1.0) print}' Arabidopsis.genotype > f_Arabidopsis.genotype
# Extract the header from the file of genotype likelihoods
head -n -1 Arabidopsis.genotype | awk '{if ($3 == NULL || $1 ~ /^@/) print}' > header_Arabidopsis.genotype
# Add footer which contains "@END_TABLE"
footer_Arabidopsis.genotype
# For large dataset: subsample a group of SNPs
sub_sample.py F_Arabidopsis.genotype -N 200000 > 200K_F_Arabidopsis.genotype
# Add the header and footer to the file of genotype likelihoods
cat header_Arabidopsis.genotype F_Arabidopsis.genotype footer_Arabidopsis.genotype > wh_wf_F_Arabidopsis.genotype
# With subsampling:
# cat header_Arabidopsis.genotype 200K_F_Arabidopsis.genotype footer_Arabidopsis.genotype > wh_wf_200K_F_Arabidopsis.genotype

# Run the relatedness command.
mapgd relatedness -i wh_wf_F_Arabidopsis.genotype -o all_Arabidopsis
