library(data.table)
counts <- fread("data/counts/ENSG00000104904-OAZ1-chr19-2270291-2270291.txt")

plot(counts$Start, counts$SRR1069690.bam)
lines(x = c(2271385, 2271442), y = c(0, 0), col = "red", lwd = 5)
lines(x = c(2271444, 2271530), y = c(0, 0), col = "red", lwd = 5)
lines(x = c(2271782, 2271954), y = c(0, 0), col = "red", lwd = 5)
lines(x = c(2272735, 2272820), y = c(0, 0), col = "red", lwd = 5)
lines(x = c(2272972, 2273490), y = c(0, 0), col = "red", lwd = 5)
