snps <- read.csv('Limb length GWAS.csv', sep = ",", header = T)
snps <- snps[, c(2,3,7)]
snps6 <- snps[snps$Pr...t.. <= 1e-6, ]
snps6 <- snps6[, c(1, 2)]
unique(snps6$CHR)

for (i in 1:7) {
  df <- snps6[snps6$CHR == paste('scaffold_',i, sep = ""), ]
  write.table(df, paste("Chr",i,"GWASe-6.txt", sep = ""), sep = "\t", quote = F, col.names = F, row.names = F)
}