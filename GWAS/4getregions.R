library(ivs)
library(dplyr)

seq <- c(1, 2, 4, 5, 6, 7)
z <- data.frame()

for (i in seq) {
  ld <- read.table(paste0("Chr",i, "GWASe-6-50Kb.geno.ld"), sep = "\t", header = T)
  ld$dist <- ld$POS2 - ld$POS1
  ld <- ld[ld$R.2 >= 0.9,]
  ld <- na.omit(ld)
  prune <- aggregate(ld$dist, by = list(ld$POS1), max)
  colnames(prune) <- c("POS1", "Distance")
  prune$CHR <- paste0("scaffold_",i)
  prune$POS2 <- prune$POS1 + prune$Distance
  pruneb <- aggregate(prune$Distance, by = list(prune$POS2), max)
  colnames(pruneb) <- c("POS2", "Distance")
  pruneb$CHR <- paste0("scaffold_",i)
  pruneb$POS1 <- pruneb$POS2 - pruneb$Distance
  pruneb <- pruneb[, c(3, 4, 1)]
  z <- rbind(z, pruneb) 
}

scaffs <- data.frame()

# Below only works for scaffs 2 and 5

for (i in seq) {
  scaff <- z[z$CHR %in% paste0('scaffold_', i),]
  sig <- scaff %>%
    mutate(Range = iv(POS1, POS2), .keep = "unused")
  regions <- iv_groups(sig$Range, abutting = T)
  regions <- as.data.frame(regions)
  x <- data.frame(Reduce(rbind, regions$regions))
  rownames(x) <- NULL
  x$CHR <- paste('scaffold_',i, sep = "")
  scaffs <- rbind(scaffs, x)
}

y <- z[z$CHR != c("scaffold_2"),]
y <- y[y$CHR != c("scaffold_5"),]
scaffs <- scaffs[, c(3, 1, 2)]
scaffs <- lapply(scaffs, unlist)
scaffs <- data.frame(scaffs)
colnames(scaffs) <- c("CHR", "POS1", "POS2")
w <- rbind(y, scaffs)

write.table(w, "GWASe-6 regions LD-pruned 0.9.txt", sep = "\t", quote = F, col.names = T, row.names = F)

# Combining LD-pruned regions with single significant SNPs, plus adding +-25Kb buffer to all
snps <- read.csv('Limb length GWAS.csv', sep = ",", header = T)
snps <- snps[,-1]

A <- read.table("GWASe-6 regions LD-pruned 0.9.txt", header = T, sep = "\t")

snps <- snps[snps$Pr...t.. <= 1e-6, ]
z <- c()
for (i in 1:nrow(snps)) {
  y <- c()
  y <- which((snps$CHR[i] == A$CHR) & (snps$POS[i] >= A$POS1) & (snps$POS[i] <= A$POS2))
  if (length(y) >= 1) {
    x <- snps[i,c(1,2)]
    z <- rbind(z,x)
  }
}

w <- c()
for (i in 1:nrow(z)) {
  y <- c()
  y <- which((z$CHR[i] == A$CHR) & (z$POS[i] >= A$POS1) & (z$POS[i] <= A$POS2))
  if (length(y) >= 1) {
    x <- z[i,c(1,2)]
    w <- rbind(w,x)
  }
}

w <- snps[, c(1,2)]
w <- na.omit(w)
k <- setdiff(w,z)

k$POS1 <- k$POS - 25000
k$POS2 <- k$POS + 25000
k <- k[, c(1,3,4)]

A$POS1 <- A$POS1 - 25000
A$POS2 <- A$POS2 + 25000

regions <- rbind(A, k)
regions <- regions[order(regions$CHR,regions$POS1 ), ]
rownames(regions) <- NULL
write.table(regions, "GWAS e-6 LD-pruned regions and single SNPs +-25kb.txt", sep = "\t", quote = F, col.names = T, row.names = F)