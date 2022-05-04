[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)

# MAPGD for plants

## Overview of an example workflow: Fastq data quality checking

MAPGD for plants.

## Installation

- __Running environment__: 
    - The workflow was constructed based on the __Linux system__.

- __Required software and versions__: 
    - [SAMtools v1.9](http://www.htslib.org/)
    - [MAPGD v0.4.38](https://github.com/LynchLab/MAPGD) 

## Input Data

The example data used here was downloaded from the *Arabidopsis* 1001 Genomes Project [JGIHeazlewood2011](https://1001genomes.org/projects/JGIHeazlewood2011/index.html).

## Submit the job to HPC cluster

```
sbatch workflow/mapgd_pipe.sh
```

## Basic steps

- Sort bam files and [filter](https://broadinstitute.github.io/picard/explain-flags.html) for mapped reads.
```
input='filelist.txt'
while IFS= read -r line
do
        echo $line
        samtools sort -o $line.sort.bam $line.bam
        samtools view -q 20 -f 3 -F 3844 -b $line.sort.bam > $line.filtered.sort.bam
        samtools index $line.filtered.sort.bam
done < $input
```

- Get header file.
```
samtools view -H Alc-0.sort.bam > Arabidopsis.header
```

- Generate mpileup file.
```
samtools mpileup -q 25 -Q 25 -B Alc-0.filtered.sort.bam Jea.filtered.sort.bam Oy-0.filtered.sort.bam Ri-0.filtered.sort.bam Sakata.filtered.sort.bam > Arabidopsis.mpileup
```

- Make a pro file of nucleotide-read quartets (counts of A, C, G, and T) from the mpileup file.
```
mapgd proview -i Arabidopsis.mpileup -H Arabidopsis.header > Arabidopsis.pro
```

- Run the allele command to estimate allele and genotype frequencies from the pro file.
```
mapgd allele -i Arabidopsis.pro -o Arabidopsis
```

- Run the genotype command to generate a file of genotype likelihoods.
```
mapgd genotype -p Arabidopsis.pro -m Arabidopsis.map > Arabidopsis.genotype
```

- Run the relatedness command.
```
mapgd relatedness -i Arabidopsis.genotype -o Arabidopsis.rel
```

## Expected results


## License
It is a free and open source software, licensed under []() (choose a license from the suggested list:  [GPLv3](https://github.com/github/choosealicense.com/blob/gh-pages/_licenses/gpl-3.0.txt), [MIT](https://github.com/github/choosealicense.com/blob/gh-pages/LICENSE.md), or [CC BY 4.0](https://github.com/github/choosealicense.com/blob/gh-pages/_licenses/cc-by-4.0.txt)).
