MAPGD output files are stored on [Figshare](https://figshare.com/s/345244bc98fa0e8a6ab1).

See [MAPGD](https://github.com/LynchLab/MAPGD) page for detailed descriptions of file types.

## filtered_Arabidopsis.map

.map files contain a list of estimated genotypic frequency obtained with the allele command. These files store test statistics for polymorphism and Hardy-Weinberg disequilibrium, as well as a small number of statistics which may prove useful for filtering variants, such as sequencing error rate and population depth of coverage.

## filtered_Arabidopsis.idx

.idx files list the name and size of all the scaffolds in a reference genome. This file can be obtained from a .bam file using the samtools view -H command and reformatting the samtools header with the 'sam2idx' command. Idx files are automatically generated when running the proview command.
	
## Arabidopsis_stats.txt	

Extracted heterozygosity scores and allele frequencies from filtered_Arabidopsis.map.

## wh_wf_F_Arabidopsis.genotype

The output of the genotype command. This stores the -log likelihood values that an individual is each of the three possible genotypes (Major Major, Major minor or minor minor) at each locus.

## header_Arabidopsis.genotype

Header file for genotype file formating.

## footer_Arabidopsis.genotype

Footer file for genotype file formating.

## Arabidopsis.rel		

The output of the relatedness command. This file stores the 7 genotypic correlation coefficients ([Ackerman et al. 2017](https://academic.oup.com/genetics/article/206/1/105/6064207)) for all pairs of individuals and some log likelihood ratio test statistics.
