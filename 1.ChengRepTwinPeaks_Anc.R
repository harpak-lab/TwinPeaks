# Cheng Rep

anc <- commandArgs(trailingOnly = TRUE)[1]
chrnum <- commandArgs(trailingOnly = TRUE)[2]

mfi <- paste0("./phasedVCF/",anc,"VCF/chr",chrnum,"_males.",anc,".INFO")
ffi <- paste0("./phasedVCF/",anc,"VCF/chr",chrnum,"_females.",anc,".INFO")

annot.genes <- read.table("annot_genes.Ensembl77.txt",sep = "\t",header = TRUE,stringsAsFactors = FALSE)

# TwinPeaksDF.filt <- read.table("TwinPeaksDF_filt_avg.txt",header = TRUE,
#                                stringsAsFactors = FALSE)
norm.counts.filt <- read.table("norm_counts_filt.Ensembl77.GTEx_v3.txt",header = TRUE,
                               stringsAsFactors = FALSE)

chr1msnps <- read.table(mfi,header = TRUE)
chr1fsnps <- read.table(ffi,header = TRUE)
chr1msnps$REF <- as.character(chr1msnps$REF)
chr1msnps$ALT <- as.character(chr1msnps$ALT)
chr1fsnps$REF <- as.character(chr1fsnps$REF)
chr1fsnps$ALT <- as.character(chr1fsnps$ALT)

chr1msnps <- chr1msnps[which(chr1msnps$POS %in% chr1fsnps$POS),]
chr1fsnps <- chr1fsnps[which(chr1fsnps$POS %in% chr1msnps$POS),]

singletons <- which(chr1msnps$AN == 1 | chr1fsnps$AN == 1)

if(length(singletons > 0)){
  chr1msnps <- chr1msnps[-singletons,]
  chr1fsnps <- chr1fsnps[-singletons,]
}

chr1msnps$AFr <- chr1msnps$AC/chr1msnps$AN
chr1fsnps$AFr <- chr1fsnps$AC/chr1fsnps$AN

chr1msnps$SE <- sqrt(chr1msnps$AN * chr1msnps$AFr * (1-chr1msnps$AFr))/chr1msnps$AN
chr1fsnps$SE <- sqrt(chr1fsnps$AN * chr1fsnps$AFr * (1-chr1fsnps$AFr))/chr1fsnps$AN

# chr1msnps$loci <- paste0("chr",chr1msnps$CHROM,":",chr1msnps$POS)
# chr1fsnps$loci <- paste0("chr",chr1fsnps$CHROM,":",chr1fsnps$POS)

FST.df.chr1 <- data.frame(NUMER = rep(NA,nrow(chr1msnps)))

getnumer <- function(pm,pf,sem,sef){
  a <- ((pm - pf)^2) - (sem^2) - (sef^2)
  return(a)
}

getdenom <- function(pm,pf,nm,nf,sem,sef){
  p <- ((pm * nm) + (pf * nf))/(nm + nf)
  b <- (4 * p * (1 - p)) - (sem^2) - (sef^2)
  return(b)
}

getnumerHud <- function(pm,pf,nm,nf){
  a <- ((pm - pf)^2)
  b <- (pm * (1 - pm))/(nm - 1)
  c <- (pf * (1 - pf))/(nf - 1)
  n <- a - b - c
  return(n)
}

getdenomHud <- function(pm,pf){
  a <- pm * (1 - pf)
  b <- pf * (1 - pm)
  d <- a + b
  return(d)
}


FST.df.chr1$NUMER <- getnumer(chr1msnps$AFr,chr1fsnps$AFr,chr1msnps$SE,chr1fsnps$SE)
FST.df.chr1$DENOM <- getdenom(chr1msnps$AFr,chr1fsnps$AFr,chr1msnps$AN,chr1fsnps$AN,chr1msnps$SE,chr1fsnps$SE)
FST.df.chr1$FST <- FST.df.chr1$NUMER/FST.df.chr1$DENOM

FST.df.chr1$NUMER_HUD <- getnumerHud(chr1msnps$AFr,chr1fsnps$AFr,chr1msnps$AN,chr1fsnps$AN)
FST.df.chr1$DENOM_HUD <- getdenomHud(chr1msnps$AFr,chr1fsnps$AFr)
FST.df.chr1$FST_HUD <- FST.df.chr1$NUMER_HUD/FST.df.chr1$DENOM_HUD

FST.df.chr1$CHROM <- chr1msnps$CHROM
FST.df.chr1$POS <- chr1msnps$POS
# FST.df.chr1$LOCI <- chr1msnps$loci

FST.df.chr1 <- FST.df.chr1[,c(7,8,1,2,3,4,5,6)]

FST.df.chr1 <- FST.df.chr1[-which(is.na(FST.df.chr1$FST)),]

chr1genes <- annot.genes$name[which(annot.genes$chr == chrnum)]

# TwinPeaksDF.filt <- norm.counts.filt[,c("GeneName","delta")]
# 
# TwinPeaksDF.chr1 <- TwinPeaksDF.filt[which(TwinPeaksDF.filt$name %in% chr1genes),]

TwinPeaksDF.rand.chr1 <- data.frame(name = norm.counts.filt$GeneName[which(norm.counts.filt$GeneName %in% chr1genes)],
                                    delta = norm.counts.filt$delta[which(norm.counts.filt$GeneName %in% chr1genes)])
TwinPeaksDF.rand.chr1$name <- as.character(TwinPeaksDF.rand.chr1$name)
TwinPeaksDF.rand.chr1$avgFST <- rep(NA,nrow(TwinPeaksDF.rand.chr1))
TwinPeaksDF.rand.chr1$maxFST <- rep(NA,nrow(TwinPeaksDF.rand.chr1))
TwinPeaksDF.rand.chr1$avgFST_HUD <- rep(NA,nrow(TwinPeaksDF.rand.chr1))
TwinPeaksDF.rand.chr1$maxFST_HUD <- rep(NA,nrow(TwinPeaksDF.rand.chr1))

GetGeneRangeRand <- function(genename){
  genestart <- (annot.genes$start[which(annot.genes$name == genename)] - 1000)
  geneend <- annot.genes$end[which(annot.genes$name == genename)]
  genechr <- annot.genes$chr[which(annot.genes$name == genename)]
  posrangel <- length(which(FST.df.chr1$POS >= genestart & 
                              FST.df.chr1$POS <= geneend &
                              FST.df.chr1$CHROM == genechr))
  if(posrangel == 0){return(0)}
  posrange <- which(FST.df.chr1$POS >= genestart & 
                      FST.df.chr1$POS <= geneend &
                      FST.df.chr1$CHROM == genechr)
  return(posrange)
}
GetGeneAvgFSTRand <- function(posrange){
  meanFST <- sum(FST.df.chr1$NUMER[posrange])/sum(FST.df.chr1$DENOM[posrange])
  return(meanFST)
}

GetGeneMaxFSTRand <- function(posrange){
  maxFST <- max(FST.df.chr1$FST[posrange])
  return(maxFST)
}

GetGeneAvgFSTHudRand <- function(posrange){
  meanFST <- mean(FST.df.chr1$FST_HUD[posrange])
  return(meanFST)
}

GetGeneMaxFSTHudRand <- function(posrange){
  maxFST <- max(FST.df.chr1$FST_HUD[posrange])
  return(maxFST)
}

for(i in 1:nrow(TwinPeaksDF.rand.chr1)){
  # cat("\r",i,"/",nrow(TwinPeaksDF.rand.chr1),sep = "")
  rangei <- GetGeneRangeRand(TwinPeaksDF.rand.chr1$name[i])
  # print(rangei)
  if(rangei == 0) next
  # print("moving on")
  TwinPeaksDF.rand.chr1$avgFST[i] <- GetGeneAvgFSTRand(rangei)
  TwinPeaksDF.rand.chr1$maxFST[i] <- GetGeneMaxFSTRand(rangei)
  TwinPeaksDF.rand.chr1$avgFST_HUD[i] <- GetGeneAvgFSTHudRand(rangei)
  TwinPeaksDF.rand.chr1$maxFST_HUD[i] <- GetGeneMaxFSTHudRand(rangei)
}

fo <- paste0("./phasedVCF/",anc,"VCF/TwinPeaksDF_chr",chrnum,".",anc,".77v3.txt")
write.table(TwinPeaksDF.rand.chr1,fo,col.names = FALSE,row.names = FALSE,
            quote = FALSE,sep = "\t")
