library(vcfR)

vcf <- read.vcfR('maf.imputed.all.vcf', verbose = T)
ext <- extract.gt(vcf)
gt <- t(ext)
gt <- as.data.frame(gt)

vec <- c("0|0" = 0, "0|1" = 1, "1|0" = 1, "1|1" = 2)
gt[] <- lapply(gt, function(x) { vec[x]})
vars <- colnames(gt)

gt$SampleID <- rownames(gt)

pheno <- read.csv("S26 Phenotypes.csv", header = T)
pca <- read.csv("maf.imputed.all.eigenvec", sep = "\t", header = T)

for (i in pheno$FID) {
  gt$relRHL.SVL.Sex[gt$SampleID == i] <- pheno$relRHL.SVL.Sex[pheno$FID == i]  
}

for (i in pca$FID) {
  gt$PC1[gt$SampleID == i] <- pca$PC1[pca$FID == i]  
}

for (i in pca$FID) {
  gt$PC2[gt$SampleID == i] <- pca$PC2[pca$FID == i]  
}

for (i in pca$FID) {
  gt$PC3[gt$SampleID == i] <- pca$PC3[pca$FID == i]  
}

for (i in pca$FID) {
  gt$PC4[gt$SampleID == i] <- pca$PC4[pca$FID == i]  
}


lm_results <- lapply(vars, function(col){
  lm_formula <- as.formula(paste("relRHL.SVL.Sex ~ PC1 + PC2 + PC3 + PC4 +", col))
  lm(lm_formula, data = gt)
})


df <- lapply(lm_results, function(x) summary(x)$coefficients)
df1 <- do.call(rbind.data.frame, df)
df2 <- df1[- grep("PC", row.names(df1)),]
df3 <- df2[- grep("Interc", row.names(df2)),]

write.csv(df3, "Limb length GWAS.csv")