library(regioneR)
geno <- read.csv("genome.csv", header = T)

A <- read.table("GWAS e-6 LD-pruned regions and single SNPs +-25kb.txt", header = T, sep = "\t")

B <- read.table("AL Tajima's D top 0.5th regions.csv", header = T, sep = ",")
B <- B[,-1]

A <- toGRanges(A)
B <- toGRanges(B)
numOverlaps(A,B)

pt <- permTest(A=A, ntimes=10000, alternative = "auto", genome = geno, randomize.function=randomizeRegions,
               evaluate.function=numOverlaps, B=B, verbose=FALSE)

plot(pt)