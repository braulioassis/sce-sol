library(ivs)
library(dplyr)

### Get significant regions for each statistic
## Alabama
# saltilassi regions

A <- read.table("SD_ALL.win200.winstep100.saltilassi.out", header = T, sep = "\t")
A <- A[A$SD_L >= quantile(A$SD_L, 0.995), ]
A$chr <- paste("scaffold_", A$chr, sep = "")
A <- A[, c(1, 2, 3)]
unique(A$chr)

for (i in 1:11) {
  las <- A[A$chr %in% paste('scaffold_',i, sep = ""),]
  if (i %in% c(7, 9, 10)){
    next
  }
  
  # Create range from windows 
  las <- las %>%
    mutate(Range = iv(start, end), .keep = "unused")
  
  # Look for overlapping windows and creating regions, list output
  regions <- iv_groups(las$Range, abutting = T)
  
  # Turning list output into dataframe
  regions <- as.data.frame(regions)
  
  x <- data.frame(Reduce(rbind, regions$regions)) #turn list into data frame
  
  rownames(x) <- NULL #remove row names 
  x$CHR <- paste('scaffold_',i, sep = "")
  
  assign(paste('scaff', i, sep = ""), x)
}

scaff8 <- A[A$chr %in% 'scaffold_8',]
scaff8$start <- scaff8$start
scaff8$end <- scaff8$end
scaff8$CHR <- scaff8$chr
scaff8 <- scaff8[, c(4, 2, 3)]

scaff11 <- A[A$chr %in% 'scaffold_11',]
scaff11$start <- scaff11$start
scaff11$end <- scaff11$end
scaff11$CHR <- scaff11$chr
scaff11 <- scaff11[, c(4, 2, 3)]

scaffs <- rbind(scaff1, scaff2, scaff3, scaff4, scaff5, scaff6, scaff8, scaff11)

scaffs <- apply(scaffs,2,as.character)
write.csv(scaffs, "AL lassi top 99.5th regions.csv")

# Tajima's D regions

A <- read.table("SD_output.100kb.20kb.TajimaD", header = T, sep = "\t")
A <- A[A$TajimaD <= quantile(A$TajimaD,0.005),]
unique(A$CHROM)

for (i in 1:22) {
  taj <- A[A$CHROM %in% paste('scaffold_',i, sep = ""),]
  if (i %in% c(7, 9)){
    next
  }
  
  # Create range from windows 
  taj <- taj %>%
    mutate(Range = iv(BIN_START, BIN_END), .keep = "unused")
  
  # Look for overlapping windows and creating regions, list output
  regions <- iv_groups(taj$Range, abutting = T)
  
  # Turning list output into dataframe
  regions <- as.data.frame(regions)
  
  x <- data.frame(Reduce(rbind, regions$regions)) #turn list into data frame
  
  rownames(x) <- NULL #remove row names 
  x$CHR <- paste('scaffold_',i, sep = "")
  
  assign(paste('scaff', i, sep = ""), x)
}

scaff10 <- A[A$CHROM %in% 'scaffold_10',]
scaff10$start <- scaff10$BIN_START
scaff10$end <- scaff10$BIN_END
scaff10$CHR <- scaff10$CHROM
scaff10 <- scaff10[, c(7, 8, 9)]

scaff19 <- A[A$CHROM %in% 'scaffold_19',]
scaff19$start <- scaff19$BIN_START
scaff19$end <- scaff19$BIN_END
scaff19$CHR <- scaff19$CHROM
scaff19 <- scaff19[, c(7, 8, 9)]

scaff20 <- A[A$CHROM %in% 'scaffold_20',]
scaff20$start <- scaff20$BIN_START
scaff20$end <- scaff20$BIN_END
scaff20$CHR <- scaff20$CHROM
scaff20 <- scaff20[, c(7, 8, 9)]

scaff21 <- A[A$CHROM %in% 'scaffold_21',]
scaff21$start <- scaff21$BIN_START
scaff21$end <- scaff21$BIN_END
scaff21$CHR <- scaff21$CHROM
scaff21 <- scaff21[, c(7, 8, 9)]

scaffs <- rbind(scaff1, scaff2, scaff3, scaff4, scaff5, scaff6, scaff8, scaff10, scaff11, scaff12, scaff13, scaff14, scaff15, scaff16, scaff17, scaff18, scaff19, scaff20, scaff21, scaff22)
scaffs <- scaffs[, c(3, 1, 2)]

scaffs <- apply(scaffs,2,as.character)
write.csv(scaffs, "AL Tajima's D top 0.5th regions.csv")

## Tennessee
# Tajima's D

A <- read.table("EE_output.100kb.20kb.TajimaD", header = T, sep = "\t")
A <- A[A$TajimaD <= quantile(A$TajimaD,0.005),]
unique(A$CHROM)

for (i in 1:22) {
  taj <- A[A$CHROM %in% paste('scaffold_',i, sep = ""),]
  if (i %in% c(9, 12, 13, 14, 17, 18, 19, 21)){
    next
  }
  
  # Create range from windows 
  taj <- taj %>%
    mutate(Range = iv(BIN_START, BIN_END), .keep = "unused")
  
  # Look for overlapping windows and creating regions, list output
  regions <- iv_groups(taj$Range, abutting = T)
  
  # Turning list output into dataframe
  regions <- as.data.frame(regions)
  
  x <- data.frame(Reduce(rbind, regions$regions)) #turn list into data frame
  
  rownames(x) <- NULL #remove row names 
  x$CHR <- paste('scaffold_',i, sep = "")
  
  assign(paste('scaff', i, sep = ""), x)
}

scaff11 <- A[A$CHROM %in% 'scaffold_11',]
scaff11$start <- scaff11$BIN_START
scaff11$end <- scaff11$BIN_END
scaff11$CHR <- scaff11$CHROM
scaff11 <- scaff11[, c(7, 8, 9)]

scaff15 <- A[A$CHROM %in% 'scaffold_15',]
scaff15$start <- scaff15$BIN_START
scaff15$end <- scaff15$BIN_END
scaff15$CHR <- scaff15$CHROM
scaff15 <- scaff15[, c(7, 8, 9)]

scaff16 <- A[A$CHROM %in% 'scaffold_16',]
scaff16$start <- scaff16$BIN_START
scaff16$end <- scaff16$BIN_END
scaff16$CHR <- scaff16$CHROM
scaff16 <- scaff16[, c(7, 8, 9)]

scaff20 <- A[A$CHROM %in% 'scaffold_20',]
scaff20$start <- scaff20$BIN_START
scaff20$end <- scaff20$BIN_END
scaff20$CHR <- scaff20$CHROM
scaff20 <- scaff20[, c(7, 8, 9)]

scaff22 <- A[A$CHROM %in% 'scaffold_22',]
scaff22$start <- scaff22$BIN_START
scaff22$end <- scaff22$BIN_END
scaff22$CHR <- scaff22$CHROM
scaff22 <- scaff22[, c(7, 8, 9)]

scaffs <- rbind(scaff1, scaff2, scaff3, scaff4, scaff5, scaff6, scaff7, scaff8, scaff10, scaff11, scaff15, scaff16, scaff20, scaff22)
scaffs <- scaffs[, c(3, 1, 2)]

scaffs <- apply(scaffs,2,as.character)
write.csv(scaffs, "TN Tajima's D top 0.5th regions.csv")

# lassi regions

A <- read.table("EE_ALL.win200.winstep100.saltilassi.out", header = T, sep = "\t")
A <- A[A$EE_L >= quantile(A$EE_L, 0.995), ]
A$chr <- paste("scaffold_", A$chr, sep = "")
A <- A[, c(1, 2, 3)]
unique(A$chr)

for (i in 1:7) {
  las <- A[A$chr %in% paste('scaffold_',i, sep = ""),]
  if (i %in% c(2, 3, 4, 5)){
    next
  }
  
  # Create range from windows 
  las <- las %>%
    mutate(Range = iv(start, end), .keep = "unused")
  
  # Look for overlapping windows and creating regions, list output
  regions <- iv_groups(las$Range, abutting = T)
  
  # Turning list output into dataframe
  regions <- as.data.frame(regions)
  
  x <- data.frame(Reduce(rbind, regions$regions)) #turn list into data frame
  
  rownames(x) <- NULL #remove row names 
  x$CHR <- paste('scaffold_',i, sep = "")
  
  assign(paste('scaff', i, sep = ""), x)
}

scaff6 <- A[A$chr %in% 'scaffold_6',]
scaff6$start <- scaff6$start
scaff6$end <- scaff6$end
scaff6$CHR <- scaff6$chr
scaff6 <- scaff6[, c(4, 2, 3)]

scaffs <- rbind(scaff1, scaff6, scaff7)

scaffs <- apply(scaffs,2,as.character)
write.csv(scaffs, "TN lassi top 99.5th regions.csv")

## Arkansas
# Tajima's D

A <- read.table("SF_output.100kb.20kb.TajimaD", header = T, sep = "\t")
A <- A[A$TajimaD <= quantile(A$TajimaD,0.005),]
unique(A$CHROM)

for (i in 1:10) {
  taj <- A[A$CHROM %in% paste('scaffold_',i, sep = ""),]
  
  # Create range from windows 
  taj <- taj %>%
    mutate(Range = iv(BIN_START, BIN_END), .keep = "unused")
  
  # Look for overlapping windows and creating regions, list output
  regions <- iv_groups(taj$Range, abutting = T)
  
  # Turning list output into dataframe
  regions <- as.data.frame(regions)
  
  x <- data.frame(Reduce(rbind, regions$regions)) #turn list into data frame
  
  rownames(x) <- NULL #remove row names 
  x$CHR <- paste('scaffold_',i, sep = "")
  
  assign(paste('scaff', i, sep = ""), x)
}

scaff9 <- A[A$CHROM %in% 'scaffold_9',]
scaff9$start <- scaff9$BIN_START
scaff9$end <- scaff9$BIN_END
scaff9$CHR <- scaff9$CHROM
scaff9 <- scaff9[, c(7, 8, 9)]

scaffs <- rbind(scaff1, scaff2, scaff3, scaff4, scaff5, scaff6, scaff7, scaff8, scaff9, scaff10)
scaffs <- scaffs[, c(3, 1, 2)]

scaffs <- apply(scaffs,2,as.character)
write.csv(scaffs, "AR Tajima's D top 0.5th regions.csv")

# lassi regions

A <- read.table("SF_ALL.win200.winstep100.saltilassi.out", header = T, sep = "\t")
A <- A[A$SF_L >= quantile(A$SF_L, 0.995), ]
A$chr <- paste("scaffold_", A$chr, sep = "")
A <- A[, c(1, 2, 3)]
unique(A$chr)

for (i in 1:6) {
  las <- A[A$chr %in% paste('scaffold_',i, sep = ""),]
  
  # Create range from windows 
  las <- las %>%
    mutate(Range = iv(start, end), .keep = "unused")
  
  # Look for overlapping windows and creating regions, list output
  regions <- iv_groups(las$Range, abutting = T)
  
  # Turning list output into dataframe
  regions <- as.data.frame(regions)
  
  x <- data.frame(Reduce(rbind, regions$regions)) #turn list into data frame
  
  rownames(x) <- NULL #remove row names 
  x$CHR <- paste('scaffold_',i, sep = "")
  
  assign(paste('scaff', i, sep = ""), x)
}

scaff4 <- A[A$chr %in% 'scaffold_4',]
scaff4$start <- scaff4$start
scaff4$end <- scaff4$end
scaff4$CHR <- scaff4$chr
scaff4 <- scaff4[, c(4, 2, 3)]

scaffs <- rbind(scaff1, scaff2, scaff3, scaff4, scaff5, scaff6)

scaffs <- apply(scaffs,2,as.character)
write.csv(scaffs, "AR lassi top 99.5th regions.csv")


### Regions significant in two of the three statistics
## Alabama
# Tajima's D and saltilassi
A <- read.csv("AL Tajima's D top 0.5th regions.csv")
A <- A[, -1]
B <- read.csv("AL lassi top 99.5th regions.csv")
B <- B[, c(4, 2, 3)]

TajLassi <- c()
for (i in 1:nrow(A)) {
  y <- c()
  y <- which(
    ((A$CHR[i] == B$CHR)&(A$start[i] >= B$start)&(A$start[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$end[i] >= B$start)&(A$end[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$start[i] <= B$start)&(A$end[i] >= B$end))
  )
  if (length(y) >= 1) {
    x <- A[i,]
    TajLassi <- rbind(TajLassi,x)
  }
}

# Tajima's D and LSBL

A <- read.csv("AL Tajima's D top 0.5th regions.csv")
A <- A[, -1]
B <- read.table("SD LSBL top0.1% LD-pruned 0.9 3+ SNPs.txt", header = T)

TajLSBL <- c()
for (i in 1:nrow(A)) {
  y <- c()
  y <- which(
    ((A$CHR[i] == B$CHR)&(A$start[i] >= B$start)&(A$start[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$end[i] >= B$start)&(A$end[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$start[i] <= B$start)&(A$end[i] >= B$end))
  )
  if (length(y) >= 1) {
    x <- A[i,]
    TajLSBL <- rbind(TajLSBL,x)
  }
}

# saltilassi and LSBL

A <- read.csv("AL lassi top 99.5th regions.csv")
A <- A[, c(4, 2, 3)]
B <- read.table("SD LSBL top0.1% LD-pruned 0.9 3+ SNPs.txt", header = T)

LassiLSBL <- c()
for (i in 1:nrow(A)) {
  y <- c()
  y <- which(
    ((A$CHR[i] == B$CHR)&(A$start[i] >= B$start)&(A$start[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$end[i] >= B$start)&(A$end[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$start[i] <= B$start)&(A$end[i] >= B$end))
  )
  if (length(y) >= 1) {
    x <- A[i,]
    LassiLSBL <- rbind(LassiLSBL,x)
  }
}

TwoOfThree <- rbind(LassiLSBL, TajLassi, TajLSBL)

unique(TwoOfThree$CHR)
seq <- c(1, 2, 3, 4, 5, 6)
scaffs <- data.frame()

for (i in seq) {
  scaff <- TwoOfThree[TwoOfThree$CHR %in% paste0('scaffold_',i),]
  # Create range from windows 
  sig <- scaff %>%
    mutate(Range = iv(start, end), .keep = "unused")
  
  # Look for overlapping windows and creating regions, list output
  regions <- iv_groups(sig$Range, abutting = T)
  
  # Turning list output into dataframe
  regions <- as.data.frame(regions)
  
  x <- data.frame(Reduce(rbind, regions$regions)) #turn list into data frame
  
  rownames(x) <- NULL #remove row names 
  x$CHR <- paste0('scaffold_',i)
  scaffs <- rbind(scaffs, x)
}

scaff8 <- data.frame(start = 36880000, end = 37180000, CHR = "scaffold_8")
scaff11 <- data.frame(start = 12200000, end = 12460000, CHR = "scaffold_11")

scaffs <- rbind(scaffs, scaff8, scaff11)

scaffs <- lapply(scaffs, unlist)
scaffs <- data.frame(scaffs)
scaffs <- scaffs[,c(3, 1, 2)]
write.table(scaffs, "AL Two out of three regions.txt", sep = "\t", quote = F, col.names = T, row.names = F)

## Tennessee
# Tajima's D and saltilassi
A <- read.csv("TN Tajima's D top 0.5th regions.csv")
A <- A[, -1]
B <- read.csv("TN lassi top 99.5th regions.csv")
B <- B[, c(4, 2, 3)]

TajLassi <- c()
for (i in 1:nrow(A)) {
  y <- c()
  y <- which(
    ((A$CHR[i] == B$CHR)&(A$start[i] >= B$start)&(A$start[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$end[i] >= B$start)&(A$end[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$start[i] <= B$start)&(A$end[i] >= B$end))
  )
  if (length(y) >= 1) {
    x <- A[i,]
    TajLassi <- rbind(TajLassi,x)
  }
}

# Tajima's D and LSBL

A <- read.csv("TN Tajima's D top 0.5th regions.csv")
A <- A[, -1]
B <- read.table("EE LSBL top0.1% LD-pruned 0.9 3+ SNPs.txt", header = T)

TajLSBL <- c()
for (i in 1:nrow(A)) {
  y <- c()
  y <- which(
    ((A$CHR[i] == B$CHR)&(A$start[i] >= B$start)&(A$start[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$end[i] >= B$start)&(A$end[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$start[i] <= B$start)&(A$end[i] >= B$end))
  )
  if (length(y) >= 1) {
    x <- A[i,]
    TajLSBL <- rbind(TajLSBL,x)
  }
}

# saltilassi and LSBL

A <- read.csv("TN lassi top 99.5th regions.csv")
A <- A[, c(4, 2, 3)]
B <- read.table("EE LSBL top0.1% LD-pruned 0.9 3+ SNPs.txt", header = T)

LassiLSBL <- c()
for (i in 1:nrow(A)) {
  y <- c()
  y <- which(
    ((A$CHR[i] == B$CHR)&(A$start[i] >= B$start)&(A$start[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$end[i] >= B$start)&(A$end[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$start[i] <= B$start)&(A$end[i] >= B$end))
  )
  if (length(y) >= 1) {
    x <- A[i,]
    LassiLSBL <- rbind(LassiLSBL,x)
  }
}

TwoOfThree <- rbind(LassiLSBL, TajLassi, TajLSBL)

unique(TwoOfThree$CHR)
seq <- c(1, 5, 6, 7)
scaffs <- data.frame()

for (i in seq) {
  scaff <- TwoOfThree[TwoOfThree$CHR %in% paste0('scaffold_',i),]
  # Create range from windows 
  sig <- scaff %>%
    mutate(Range = iv(start, end), .keep = "unused")
  
  # Look for overlapping windows and creating regions, list output
  regions <- iv_groups(sig$Range, abutting = T)
  
  # Turning list output into dataframe
  regions <- as.data.frame(regions)
  
  x <- data.frame(Reduce(rbind, regions$regions)) #turn list into data frame
  
  rownames(x) <- NULL #remove row names 
  x$CHR <- paste0('scaffold_',i)
  scaffs <- rbind(scaffs, x)
}

scaff2 <- data.frame(start = 299040000, end = 299300000, CHR = "scaffold_2")
scaff3 <- data.frame(start = 216560000, end = 216800000, CHR = "scaffold_3")
scaff8 <- data.frame(start = 24980000, end = 25160000, CHR = "scaffold_8")

scaffs <- rbind(scaffs, scaff2, scaff3, scaff8)

scaffs <- lapply(scaffs, unlist)
scaffs <- data.frame(scaffs)
scaffs <- scaffs[,c(3, 1, 2)]
scaffs <- scaffs[order(scaffs$CHR, scaffs$start), ]
write.table(scaffs, "TN Two out of three regions.txt", sep = "\t", quote = F, col.names = T, row.names = F)

## Arkansas
# Tajima's D and saltilassi
A <- read.csv("AR Tajima's D top 0.5th regions.csv")
A <- A[, -1]
B <- read.csv("AR lassi top 99.5th regions.csv")
B <- B[, c(4, 2, 3)]

TajLassi <- c()
for (i in 1:nrow(A)) {
  y <- c()
  y <- which(
    ((A$CHR[i] == B$CHR)&(A$start[i] >= B$start)&(A$start[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$end[i] >= B$start)&(A$end[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$start[i] <= B$start)&(A$end[i] >= B$end))
  )
  if (length(y) >= 1) {
    x <- A[i,]
    TajLassi <- rbind(TajLassi,x)
  }
}

# Tajima's D and LSBL

A <- read.csv("AR Tajima's D top 0.5th regions.csv")
A <- A[, -1]
B <- read.table("SF LSBL top0.1% LD-pruned 0.9 3+ SNPs.txt", header = T)

TajLSBL <- c()
for (i in 1:nrow(A)) {
  y <- c()
  y <- which(
    ((A$CHR[i] == B$CHR)&(A$start[i] >= B$start)&(A$start[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$end[i] >= B$start)&(A$end[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$start[i] <= B$start)&(A$end[i] >= B$end))
  )
  if (length(y) >= 1) {
    x <- A[i,]
    TajLSBL <- rbind(TajLSBL,x)
  }
}

# saltilassi and LSBL

A <- read.csv("AR lassi top 99.5th regions.csv")
A <- A[, c(4, 2, 3)]
B <- read.table("SF LSBL top0.1% LD-pruned 0.9 3+ SNPs.txt", header = T)

LassiLSBL <- c()
for (i in 1:nrow(A)) {
  y <- c()
  y <- which(
    ((A$CHR[i] == B$CHR)&(A$start[i] >= B$start)&(A$start[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$end[i] >= B$start)&(A$end[i] <= B$end))|
      ((A$CHR[i] == B$CHR)&(A$start[i] <= B$start)&(A$end[i] >= B$end))
  )
  if (length(y) >= 1) {
    x <- A[i,]
    LassiLSBL <- rbind(LassiLSBL,x)
  }
}

TwoOfThree <- rbind(LassiLSBL, TajLassi, TajLSBL)

unique(TwoOfThree$CHR)
seq <- c(1, 2, 3, 5, 6)
scaffs <- data.frame()

for (i in seq) {
  scaff <- TwoOfThree[TwoOfThree$CHR %in% paste0('scaffold_',i),]
  # Create range from windows 
  sig <- scaff %>%
    mutate(Range = iv(start, end), .keep = "unused")
  
  # Look for overlapping windows and creating regions, list output
  regions <- iv_groups(sig$Range, abutting = T)
  
  # Turning list output into dataframe
  regions <- as.data.frame(regions)
  
  x <- data.frame(Reduce(rbind, regions$regions)) #turn list into data frame
  
  rownames(x) <- NULL #remove row names 
  x$CHR <- paste0('scaffold_',i)
  scaffs <- rbind(scaffs, x)
}

scaff4 <- data.frame(start = 78260000, end = 78420000, CHR = "scaffold_4")
scaff8 <- data.frame(start = 9360000, end = 9480000, CHR = "scaffold_8")

scaffs <- rbind(scaffs, scaff4, scaff8)

scaffs <- lapply(scaffs, unlist)
scaffs <- data.frame(scaffs)
scaffs <- scaffs[,c(3, 1, 2)]
scaffs <- scaffs[order(scaffs$CHR, scaffs$start), ]
write.table(scaffs, "AR Two out of three regions.txt", sep = "\t", quote = F, col.names = T, row.names = F)
