# Alabama
LSBL <- read.table("EESF.weir.fst", sep = "\t", header = T)
EESD <- read.table("EESD.weir.fst", sep = "\t", header = T)
SDSF <- read.table("SDSF.weir.fst", sep = "\t", header = T)
LSBL$WeirEESF <- LSBL$WEIR_AND_COCKERHAM_FST
LSBL$WeirEESD <- EESD$WEIR_AND_COCKERHAM_FST
LSBL$WeirSDSF <- SDSF$WEIR_AND_COCKERHAM_FST

LSBL <- LSBL[,-3]

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

LSBL[is.nan(LSBL)] <- 0

LSBL$SD.LSBL <- (LSBL$WeirEESD + LSBL$WeirSDSF - LSBL$WeirEESF)/2

write.csv(LSBL,"LSBL.SD.csv")

# Arkansas
LSBL <- read.table("EESD.weir.fst", sep = "\t", header = T)
EESF <- read.table("EESF.weir.fst", sep = "\t", header = T)
SDSF <- read.table("SDSF.weir.fst", sep = "\t", header = T)
LSBL$WeirEESD <- LSBL$WEIR_AND_COCKERHAM_FST
LSBL$WeirEESF <- EESF$WEIR_AND_COCKERHAM_FST
LSBL$WeirSDSF <- SDSF$WEIR_AND_COCKERHAM_FST

LSBL <- LSBL[,-3]

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

LSBL[is.nan(LSBL)] <- 0

LSBL$SF.LSBL <- (LSBL$WeirEESF + LSBL$WeirSDSF - LSBL$WeirEESD)/2

write.csv(LSBL,"LSBL.SF.csv")

# Tennessee
LSBL <- read.table("SDSF.weir.fst", sep = "\t", header = T)
EESF <- read.table("EESF.weir.fst", sep = "\t", header = T)
EESD <- read.table("EESD.weir.fst", sep = "\t", header = T)
LSBL$WeirSDSF <- LSBL$WEIR_AND_COCKERHAM_FST
LSBL$WeirEESF <- EESF$WEIR_AND_COCKERHAM_FST
LSBL$WeirSDSF <- SDSF$WEIR_AND_COCKERHAM_FST

LSBL <- LSBL[,-3]

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

LSBL[is.nan(LSBL)] <- 0

LSBL$EE.LSBL <- (LSBL$WeirEESF + LSBL$WeirEESD - LSBL$WeirSDSF)/2

write.csv(LSBL,"LSBL.EE.csv")