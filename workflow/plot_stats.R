library("dplyr")
library(Hmisc)

Data <- read.table("Arabidopsis_stats.txt", head = T)
Data$Mn_FREQ <- with(Data, 1 - MJ_FREQ)
Chrom1 <- Data[which( Data$SCFNAME ==1 ),]

Chrom1_windows <- Chrom1 %>%                        # 100kb window
  mutate(ranges = cut(POS,
                      seq(21, 30417973, 100000))) %>% 
  group_by(ranges) %>% 
  dplyr::summarize(het = mean(HETERO), freq = mean(Mn_FREQ)) %>% 
  as.data.frame()

par(mfrow=c(2,1))
plot(x=Chrom1_windows$ranges, y=Chrom1_windows$het, xlab="", ylab="Heterozygosity",xaxt="n")
plot(x=Chrom1_windows$ranges, y=Chrom1_windows$freq, xlab="Chrom 1", ylab="Minor allele frequency",xaxt="n")
axis(side=1,at=c(1,150,300),labels=c("1","15M","30M"))