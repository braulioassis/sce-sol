#Background set of genes

annot <- read.csv("SceUnd1.0_top24.annotations.txt", sep = "\t", header =T)
annot <- annot[annot$Product != "hypothetical protein",]
annot$Name <- gsub("_\\d", "", annot$Name)
vec <- unique(annot$Name)
write(vec, "background.txt")

# For regions significant in two out of three statistics
annot <- read.csv("SceUnd1.0_top24.annotations.txt", sep = "\t", header =T)
annot <- annot[annot$Product != "hypothetical protein",]
annot$Name <- gsub("_\\d", "", annot$Name)

regions <- read.table("AL Two out of three regions.txt", header = T)

vec <- c()
for (i in 1:nrow(annot)) {
  y <- c()
  y <- which(
    ((annot$Contig[i] == regions$CHR)&(annot$Start[i] >= regions$start-25000)&(annot$Start[i] <= regions$end+25000))
    | ((annot$Contig[i] == regions$CHR)&(annot$Stop[i] >= regions$start-25000)&(annot$Stop[i] <= regions$end+25000))
    | ((annot$Contig[i] == regions$CHR)&(annot$Start[i] >= regions$start)&(annot$Stop[i] <= regions$end))
  )
  if (length(y) >= 1) {
    x <- annot[i,7]
    vec <- rbind(vec,x)
  }
}

vec <- gsub("_\\d", "", vec)
vec <- unique(vec)
write(vec, "AL genes for two out of three regions +- 25kb.txt")

# Tajima's D significant genes

annot <- read.csv("SceUnd1.0_top24.annotations.txt", sep = "\t", header =T)
annot <- annot[annot$Product != "hypothetical protein",]
annot$Name <- gsub("_\\d", "", annot$Name)

Tajima <- read.table("SD_output.100kb.20kb.TajimaD", header = T, sep = "\t")
Tajima <- Tajima[Tajima$TajimaD <= quantile(Tajima$TajimaD,0.005), ]

vec <- c()
for (i in 1:nrow(annot)) {
  y <- c()
  y <- which(
    ((annot$Contig[i] == Tajima$CHROM)&(annot$Start[i] >= Tajima$BIN_START-25000)&(annot$Start[i] <= Tajima$BIN_END+25000))
    | ((annot$Contig[i] == Tajima$CHROM)&(annot$Stop[i] >= Tajima$BIN_START-25000)&(annot$Stop[i] <= Tajima$BIN_END+25000))
    | ((annot$Contig[i] == Tajima$CHROM)&(annot$Start[i] >= Tajima$BIN_START)&(annot$Stop[i] <= Tajima$BIN_END))
  )
  if (length(y) >= 1) {
    x <- annot[i,7]
    vec <- rbind(vec,x)
  }
}

vec <- gsub("_\\d", "", vec)
vec <- unique(vec)
write(vec, "Tajima's D 100k+-25kb 0.5th percentile any gene intersecting window.txt")

# saltilassi significant genes

annot <- read.csv("SceUnd1.0_top24.annotations.txt", sep = "\t", header =T)
annot <- annot[annot$Product != "hypothetical protein",]
annot$Name <- gsub("_\\d", "", annot$Name)

lassi <- read.table("SD_ALL.win200.winstep100.saltilassi.out", header = T, sep = "\t")
lassi <- lassi[lassi$SD_L >= quantile(lassi$SD_L, 0.995), ]
lassi$chr <- paste0("scaffold_", lassi$chr)

vec <- c()
for (i in 1:nrow(annot)) {
  y <- c()
  y <- which(((annot$Contig[i] == lassi$chr)&(annot$Start[i] >= lassi$start-25000)&(annot$Start[i] <= lassi$end+25000))
             |((annot$Contig[i] == lassi$chr)&(annot$Stop[i] >= lassi$start-25000)&(annot$Stop[i] <= lassi$end+25000))
             |((annot$Contig[i] == lassi$chr)&(annot$Start[i] >= lassi$start)&(annot$Stop[i] <= lassi$end))
  )
  if (length(y) >= 1) {
    x <- annot[i,7]
    vec <- rbind(vec,x)
  }
}

vec <- gsub("_\\d", "", vec)
vec <- unique(vec)
write(vec, "lassi +-25kb  99.5th percentile any gene intersecting window.txt")

# Repeat for other populations