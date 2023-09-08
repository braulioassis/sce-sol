## For a pair populations, find all positions that are not fixed in both for the same allele, and
## all positions that a fixed allele in one population is the major allele in the other population
library(tidyr)

# Tennessee/Alabama (for Fst)

df <- read.table("EE.frq", header = T, sep = "\t")
SD <- read.table("SD.frq", header = T, sep = "\t")

df <- separate(df, EE.a, into = c("EE.AL.a", "EE.AL.a.frq"), sep = ":", remove = FALSE)
df$EE.AL.a.frq <- gsub("\\s+", "", df$EE.AL.a.frq)

df <- separate(df, EE.b, into = c("EE.AL.b", "EE.AL.b.frq"), sep = ":", remove = FALSE)
df$EE.AL.b.frq <- gsub("\\s+", "", df$EE.AL.b.frq)

df$SD.a <- SD$SD.a
df$SD.b <- SD$SD.b
SD <- c()

df <- separate(df, SD.a, into = c("SD.AL.a", "SD.AL.a.frq"), sep = ":", remove = FALSE)
df$SD.AL.a.frq <- gsub("\\s+", "", df$SD.AL.a.frq)

df <- separate(df, SD.b, into = c("SD.AL.b", "SD.AL.b.frq"), sep = ":", remove = FALSE)
df$SD.AL.b.frq <- gsub("\\s+", "", df$SD.AL.b.frq)

# Iterate over each row
for (i in 1:nrow(df)) {
  # Check if EE.AL.a is different from SD.AL.a
  if (df[i, "EE.AL.a"] != df[i, "SD.AL.a"]) {
    # Swap values between SD.AL.a and SD.AL.b
    temp5 <- df[i, "SD.AL.a"]
    df[i, "SD.AL.a"] <- df[i, "SD.AL.b"]
    df[i, "SD.AL.b"] <- temp5
    
    # Swap values between SD.AL.a.frq and SD.AL.b.frq
    temp6 <- df[i, "SD.AL.a.frq"]
    df[i, "SD.AL.a.frq"] <- df[i, "SD.AL.b.frq"]
    df[i, "SD.AL.b.frq"] <- temp6
  }
}

write.table(nrow(which(df$EE.AL.a != df$SD.AL.a)), "diffallelesEESD.txt")

# Keep only positions that are not fixed in both populations
df <- subset(df, !(EE.AL.a.frq == 1 & SD.AL.a.frq == 1) & !(EE.AL.a.frq == 0 & SD.AL.a.frq == 0))
# Keep only positions that when one allele is fixed in one population, that allele is the major allele in the other population
df <- subset(df, !(EE.AL.a.frq == 0 & SD.AL.a.frq < 0.5) & !(EE.AL.a.frq == 1 & SD.AL.a.frq > 0.5))

write.table(df[, c(1,2)], "EE-SD-positionstokeep.txt", sep = "\t", quote = F, col.names = F, row.names = F)

# Tennessee/Arkansas (for Fst)

df <- read.table("EE.frq", header = T, sep = "\t")
SF <- read.table("SF.frq", header = T, sep = "\t")

df <- separate(df, EE.a, into = c("EE.AL.a", "EE.AL.a.frq"), sep = ":", remove = FALSE)
df$EE.AL.a.frq <- gsub("\\s+", "", df$EE.AL.a.frq)

df <- separate(df, EE.b, into = c("EE.AL.b", "EE.AL.b.frq"), sep = ":", remove = FALSE)
df$EE.AL.b.frq <- gsub("\\s+", "", df$EE.AL.b.frq)

df$SF.a <- SF$SF.a
df$SF.b <- SF$SF.b
SF <- c()

df <- separate(df, SF.a, into = c("SF.AL.a", "SF.AL.a.frq"), sep = ":", remove = FALSE)
df$SF.AL.a.frq <- gsub("\\s+", "", df$SF.AL.a.frq)

df <- separate(df, SF.b, into = c("SF.AL.b", "SF.AL.b.frq"), sep = ":", remove = FALSE)
df$SF.AL.b.frq <- gsub("\\s+", "", df$SF.AL.b.frq)

# Iterate over each row
for (i in 1:nrow(df)) {
  # Check if EE.AL.a is different from SF.AL.a
  if (df[i, "EE.AL.a"] != df[i, "SF.AL.a"]) {
    # Swap values between SF.AL.a and SF.AL.b
    temp5 <- df[i, "SF.AL.a"]
    df[i, "SF.AL.a"] <- df[i, "SF.AL.b"]
    df[i, "SF.AL.b"] <- temp5
    
    # Swap values between SF.AL.a.frq and SF.AL.b.frq
    temp6 <- df[i, "SF.AL.a.frq"]
    df[i, "SF.AL.a.frq"] <- df[i, "SF.AL.b.frq"]
    df[i, "SF.AL.b.frq"] <- temp6
  }
}

write.table(nrow(which(df$EE.AL.a != df$SF.AL.a)), "diffallelesEESF.txt")

# Keep only positions that are not fixed in both populations
df <- subset(df, !(EE.AL.a.frq == 1 & SF.AL.a.frq == 1) & !(EE.AL.a.frq == 0 & SF.AL.a.frq == 0))
# Keep only positions that when one allele is fixed in one population, that allele is the major allele in the other population
df <- subset(df, !(EE.AL.a.frq == 0 & SF.AL.a.frq < 0.5) & !(EE.AL.a.frq == 1 & SF.AL.a.frq > 0.5))

write.table(df[, c(1,2)], "EE-SF-positionstokeep.txt", sep = "\t", quote = F, col.names = F, row.names = F)

# Arkansas/Alabama (for Fst)

df <- read.table("SF.frq", header = T, sep = "\t")
SD <- read.table("SD.frq", header = T, sep = "\t")

df <- separate(df, SF.a, into = c("SF.AL.a", "SF.AL.a.frq"), sep = ":", remove = FALSE)
df$SF.AL.a.frq <- gsub("\\s+", "", df$SF.AL.a.frq)

df <- separate(df, SF.b, into = c("SF.AL.b", "SF.AL.b.frq"), sep = ":", remove = FALSE)
df$SF.AL.b.frq <- gsub("\\s+", "", df$SF.AL.b.frq)

df$SD.a <- SD$SD.a
df$SD.b <- SD$SD.b
SD <- c()

df <- separate(df, SD.a, into = c("SD.AL.a", "SD.AL.a.frq"), sep = ":", remove = FALSE)
df$SD.AL.a.frq <- gsub("\\s+", "", df$SD.AL.a.frq)

df <- separate(df, SD.b, into = c("SD.AL.b", "SD.AL.b.frq"), sep = ":", remove = FALSE)
df$SD.AL.b.frq <- gsub("\\s+", "", df$SD.AL.b.frq)

# Iterate over each row
for (i in 1:nrow(df)) {
  # Check if SF.AL.a is different from SD.AL.a
  if (df[i, "SF.AL.a"] != df[i, "SD.AL.a"]) {
    # Swap values between SD.AL.a and SD.AL.b
    temp5 <- df[i, "SD.AL.a"]
    df[i, "SD.AL.a"] <- df[i, "SD.AL.b"]
    df[i, "SD.AL.b"] <- temp5
    
    # Swap values between SD.AL.a.frq and SD.AL.b.frq
    temp6 <- df[i, "SD.AL.a.frq"]
    df[i, "SD.AL.a.frq"] <- df[i, "SD.AL.b.frq"]
    df[i, "SD.AL.b.frq"] <- temp6
  }
}

write.table(nrow(which(df$SF.AL.a != df$SD.AL.a)), "diffallelesSFSD.txt")

# Keep only positions that are not fixed in both populations
df <- subset(df, !(SF.AL.a.frq == 1 & SD.AL.a.frq == 1) & !(SF.AL.a.frq == 0 & SD.AL.a.frq == 0))
# Keep only positions that when one allele is fixed in one population, that allele is the major allele in the other population
df <- subset(df, !(SF.AL.a.frq == 0 & SD.AL.a.frq < 0.5) & !(SF.AL.a.frq == 1 & SD.AL.a.frq > 0.5))

write.table(df[, c(1,2)], "SF-SD-positionstokeep.txt", sep = "\t", quote = F, col.names = F, row.names = F)


## For the three populations, remove all positions where one allele is fixed in two populations,
## and that same allele is the minor allele in the third population

df <- read.table("SD.frq", header = T, sep = "\t")
SF <- read.table("SF.frq", header = T, sep = "\t")
EE <- read.table("EE.frq", header = T, sep = "\t")

df <- separate(df, SD.a, into = c("SD.AL.a", "SD.AL.a.frq"), sep = ":", remove = FALSE)
df$SD.AL.a.frq <- gsub("\\s+", "", df$SD.AL.a.frq)

df <- separate(df, SD.b, into = c("SD.AL.b", "SD.AL.b.frq"), sep = ":", remove = FALSE)
df$SD.AL.b.frq <- gsub("\\s+", "", df$SD.AL.b.frq)

df$SF.a <- SF$SF.a
df$SF.b <- SF$SF.b

df <- separate(df, SF.a, into = c("SF.AL.a", "SF.AL.a.frq"), sep = ":", remove = FALSE)
df$SF.AL.a.frq <- gsub("\\s+", "", df$SF.AL.a.frq)

df <- separate(df, SF.b, into = c("SF.AL.b", "SF.AL.b.frq"), sep = ":", remove = FALSE)
df$SF.AL.b.frq <- gsub("\\s+", "", df$SF.AL.b.frq)

df$EE.a <- EE$EE.a
df$EE.b <- EE$EE.b

df <- separate(df, EE.a, into = c("EE.AL.a", "EE.AL.a.frq"), sep = ":", remove = FALSE)
df$EE.AL.a.frq <- gsub("\\s+", "", df$EE.AL.a.frq)

df <- separate(df, EE.b, into = c("EE.AL.b", "EE.AL.b.frq"), sep = ":", remove = FALSE)
df$EE.AL.b.frq <- gsub("\\s+", "", df$EE.AL.b.frq)

SF <- c()
EE <- c()

# Iterate over each row
for (i in 1:nrow(df)) {
  # Check if SD.AL.a is different from SF.AL.a
  if (df[i, "SD.AL.a"] != df[i, "SF.AL.a"]) {
    # Swap values between SF.AL.a and SF.AL.b
    temp5 <- df[i, "SF.AL.a"]
    df[i, "SF.AL.a"] <- df[i, "SF.AL.b"]
    df[i, "SF.AL.b"] <- temp5
    
    # Swap values between SF.AL.a.frq and SF.AL.b.frq
    temp6 <- df[i, "SF.AL.a.frq"]
    df[i, "SF.AL.a.frq"] <- df[i, "SF.AL.b.frq"]
    df[i, "SF.AL.b.frq"] <- temp6
  }
}

for (i in 1:nrow(df)) {
  # Check if SF.AL.a is different from EE.AL.a
  if (df[i, "SF.AL.a"] != df[i, "EE.AL.a"]) {
    # Swap values between EE.AL.a and EE.AL.b
    temp5 <- df[i, "EE.AL.a"]
    df[i, "EE.AL.a"] <- df[i, "EE.AL.b"]
    df[i, "EE.AL.b"] <- temp5
    
    # Swap values between EE.AL.a.frq and EE.AL.b.frq
    temp6 <- df[i, "EE.AL.a.frq"]
    df[i, "EE.AL.a.frq"] <- df[i, "EE.AL.b.frq"]
    df[i, "EE.AL.b.frq"] <- temp6
  }
}

write.table(nrow(which(df$SD.AL.a != df$SF.AL.a)), "diffallelesSDSF.txt")
write.table(nrow(which(df$SD.AL.a != df$EE.AL.a)), "diffallelesSDEE.txt")
write.table(nrow(which(df$SF.AL.a != df$EE.AL.a)), "diffallelesSFEE.txt")

df <- subset(df, !(SD.AL.a.frq == 1 & SF.AL.a.frq == 1 & EE.AL.a.frq > 0.5) & !(SD.AL.a.frq == 0 & SF.AL.a.frq == 0 & EE.AL.a.frq < 0.5))
df <- subset(df, !(SD.AL.a.frq == 1 & EE.AL.a.frq == 1 & SF.AL.a.frq > 0.5) & !(SD.AL.a.frq == 0 & EE.AL.a.frq == 0 & SF.AL.a.frq < 0.5))
df <- subset(df, !(EE.AL.a.frq == 1 & SF.AL.a.frq == 1 & SD.AL.a.frq > 0.5) & !(EE.AL.a.frq == 0 & SF.AL.a.frq == 0 & SD.AL.a.frq < 0.5))

write.table(df[, c(1,2)], "SD-SF-EE-positionstokeep.txt", sep = "\t", quote = F, col.names = F, row.names = F)