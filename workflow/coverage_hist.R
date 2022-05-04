Coverage <- read.table("Arabidopsis_coverage.txt", head = T)
hist(Coverage$COVERAG, xlim = c(0, 400), breaks = 100000, xlab = "Coverage", main = "")
