### LSBL
## Alabama
# Get top LSBL positions
LSBL <- read.csv("LSBL.SD.csv")
LSBL <- LSBL[LSBL$SD.LSBL >= quantile(LSBL$SD.LSBL, 0.999), ]
LSBL <- LSBL[, c(2,3)]

for (i in 1:24) {
  df <- LSBL[LSBL$CHROM == paste('scaffold_',i, sep = ""),]
  write.table(df, paste("SDchr",i,"LSBLtopSNPs.txt", sep = ""), sep = "\t", quote = F, col.names = F, row.names = F)
}

## Arkansas
# Get top LSBL positions
LSBL <- read.csv("LSBL.SF.csv")
LSBL <- LSBL[LSBL$SF.LSBL >= quantile(LSBL$SF.LSBL, 0.999), ]
LSBL <- LSBL[, c(2,3)]

for (i in 1:24) {
  df <- LSBL[LSBL$CHROM == paste('scaffold_',i, sep = ""),]
  write.table(df, paste("SFchr",i,"LSBLtopSNPs.txt", sep = ""), sep = "\t", quote = F, col.names = F, row.names = F)
}

## Tennessee
# Get top LSBL positions
LSBL <- read.csv("LSBL.EE.csv")
LSBL <- LSBL[LSBL$EE.LSBL >= quantile(LSBL$EE.LSBL, 0.999), ]
LSBL <- LSBL[, c(2,3)]

for (i in 1:24) {
  df <- LSBL[LSBL$CHROM == paste('scaffold_',i, sep = ""),]
  write.table(df, paste("EEchr",i,"LSBLtopSNPs.txt", sep = ""), sep = "\t", quote = F, col.names = F, row.names = F)
}