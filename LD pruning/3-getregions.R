library(ivs)
library(dplyr)
### Get regions with 3+ SNPs in the top 0.1% of LSBL distribution
###and with r^2 > 0.9 and less than 25Kb apart

## Alabama

z <- data.frame()
seq <- c(1:12)
for (i in seq) {
  ld <- read.table(paste0("SDchr",i, "-LSBLtop0.1-50kb.geno.ld"), sep = "\t", header = T)
  pos <- read.table(paste0("SDchr",i,"LSBLtopSNPs.txt"), sep = "\t")
  ld <- ld[ld$POS2 %in% pos$V2,]
  if (nrow(ld) == 0) {
    print(paste0("No linked SNPs in Chr ", i))
  } else {
    ld$dist <- ld$POS2 - ld$POS1
    ld <- ld[ld$R.2 >= 0.9,]
    ld <- na.omit(ld)
    prune <- aggregate(ld$dist, by = list(ld$POS1), max)
    colnames(prune) <- c("POS1", "Distance")
    prune$CHR <- paste0("scaffold_",i)
    prune$POS2 <- prune$POS1 + prune$Distance
    prune <- prune[, c(3, 1, 4)]
    z <- rbind(z, prune) 
  }
}

scaffs <- data.frame()

unique(z$CHR)
for (i in 1:12) {
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
unique(scaffs$CHR)

# Bind scaff 12 manually if needed
scaff12 <- data.frame(start = 3806591, end = 3807049, CHR = "scaffold_12")
scaffs <- rbind(scaffs,scaff12)

scaffs <- scaffs[,c(3, 1, 2)]
scaffs <- lapply(scaffs, unlist)
scaffs <- data.frame(scaffs)
write.table(scaffs, "SD LSBL top0.1% LD-pruned 0.9.txt", sep = "\t", quote = F, col.names = T, row.names = F)

## Arkansas

z <- data.frame()
seq <- c(1:12, 15, 16, 19)
for (i in seq) {
  ld <- read.table(paste0("SFchr",i, "-LSBLtop0.1-50kb.geno.ld"), sep = "\t", header = T)
  pos <- read.table(paste0("SFchr",i,"LSBLtopSNPs.txt"), sep = "\t")
  ld <- ld[ld$POS2 %in% pos$V2,]
  if (nrow(ld) == 0) {
    print(paste0("No linked SNPs in Chr ", i))
  } else {
    ld$dist <- ld$POS2 - ld$POS1
    ld <- ld[ld$R.2 >= 0.9,]
    ld <- na.omit(ld)
    prune <- aggregate(ld$dist, by = list(ld$POS1), max)
    colnames(prune) <- c("POS1", "Distance")
    prune$CHR <- paste0("scaffold_",i)
    prune$POS2 <- prune$POS1 + prune$Distance
    prune <- prune[, c(3, 1, 4)]
    z <- rbind(z, prune) 
  }
}

scaffs <- data.frame()

unique(z$CHR)
for (i in 1:12) {
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
unique(scaffs$CHR)

# Bind scaff 12 manually if needed
scaff12 <- data.frame(start = 1561123, end = 1584680, CHR = "scaffold_12")
scaffs <- rbind(scaffs,scaff12)

scaffs <- scaffs[,c(3, 1, 2)]
scaffs <- lapply(scaffs, unlist)
scaffs <- data.frame(scaffs)
write.table(scaffs, "SF LSBL top0.1% LD-pruned 0.9.txt", sep = "\t", quote = F, col.names = T, row.names = F)

## Tennessee

z <- data.frame()
seq <- c(1:12, 14:17, 19)
for (i in seq) {
  ld <- read.table(paste0("EEchr",i, "-LSBLtop0.1-50kb.geno.ld"), sep = "\t", header = T)
  pos <- read.table(paste0("EEchr",i,"LSBLtopSNPs.txt"), sep = "\t")
  ld <- ld[ld$POS2 %in% pos$V2,]
  if (nrow(ld) == 0) {
    print(paste0("No linked SNPs in Chr ", i))
  } else {
    ld$dist <- ld$POS2 - ld$POS1
    ld <- ld[ld$R.2 >= 0.9,]
    ld <- na.omit(ld)
    prune <- aggregate(ld$dist, by = list(ld$POS1), max)
    colnames(prune) <- c("POS1", "Distance")
    prune$CHR <- paste0("scaffold_",i)
    prune$POS2 <- prune$POS1 + prune$Distance
    prune <- prune[, c(3, 1, 4)]
    z <- rbind(z, prune) 
  }
}

scaffs <- data.frame()

unique(z$CHR)
for (i in 1:11) {
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
unique(scaffs$CHR)

scaffs <- scaffs[,c(3, 1, 2)]
scaffs <- lapply(scaffs, unlist)
scaffs <- data.frame(scaffs)
write.table(scaffs, "EE LSBL top0.1% LD-pruned 0.9.txt", sep = "\t", quote = F, col.names = T, row.names = F)