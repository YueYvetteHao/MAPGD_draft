library("dplyr")
library(Hmisc)

Hetero <- read.table("Arabidopsis_het.txt", head = T)
Chrom1 <- Hetero[which( Hetero$SCFNAME ==1 ),]

Chrom1_windows <- Chrom1 %>%                        # 100kb window
  mutate(ranges = cut(POS,
                      seq(21, 30417973, 100000))) %>% 
  group_by(ranges) %>% 
  dplyr::summarize(mean = mean(HETERO)) %>% 
  as.data.frame()

plot(x=Chrom1_windows$ranges, y=Chrom1_windows$mean, xlab="Chrom 1", ylab="Heterozygosity",xaxt="n")
axis(side=1,at=c(1,150,300),labels=c("1","15M","30M"))