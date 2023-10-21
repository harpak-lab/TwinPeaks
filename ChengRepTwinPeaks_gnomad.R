# ChengRep Twin Peaks gnomAD v3 + GTEx v8 + GRCh38

setwd("/scratch/08312/mjm8356/ChengRep/gnomAD/")

anc <- commandArgs(trailingOnly = TRUE)[1]
samptype <- commandArgs(trailingOnly = TRUE)[2]
chrnum <- commandArgs(trailingOnly = TRUE)[3]

exprdata <- read.table("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",
                       header = FALSE,fill = TRUE,sep = "\t")

sampleAtt <- read.table("GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                        header = TRUE,fill = TRUE,sep = "\t")

samplePhen <- read.table("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
                         header = TRUE,fill = TRUE,sep = "\t")

# Tissue Types: "Blood, Brain, Adipose Tissue, Muscle, Blood Vessel, Heart, Ovary,
# Uterus, Vagina, Breast, Skin, Salivary Gland, Adrenal Gland, Thyroid, Lung, Spleen,
# Pancreas, Esophagus, Stomach, Colon, Small Intestine, Prostate, Testis, Nerve,
# Liver, Pituitary, Kidney, "

if(samptype == "Gonads"){
  sampleAtt <- sampleAtt[which(sampleAtt$SMTS == "Ovary" | sampleAtt$SMTS == "Testis"),]
}else{
  sampleAtt <- sampleAtt[which(sampleAtt$SMTS == samptype),]
}
msubs <- samplePhen$SUBJID[which(samplePhen$SEX == 1)]
fsubs <- samplePhen$SUBJID[which(samplePhen$SEX == 2)]
msamps <- sampleAtt$SAMPID[grepl(paste(msubs,collapse = "|"),sampleAtt$SAMPID)]
fsamps <- sampleAtt$SAMPID[grepl(paste(fsubs,collapse = "|"),sampleAtt$SAMPID)]

exprdata <- exprdata[-1,]
names(exprdata) <- exprdata[1,]
exprdata <- exprdata[-1,]

DeltaTab <- exprdata[,c(1,which(names(exprdata) %in% msamps),which(names(exprdata) %in% fsamps))]
DeltaTab[,2:ncol(DeltaTab)] <- lapply(DeltaTab[,2:ncol(DeltaTab)],as.numeric)
DeltaTab$mexpr <- rowMeans(DeltaTab[1:nrow(DeltaTab),which(names(DeltaTab) %in% msamps)],na.rm = TRUE)
DeltaTab$fexpr <- rowMeans(DeltaTab[1:nrow(DeltaTab),which(names(DeltaTab) %in% fsamps)],na.rm = TRUE)
DeltaTab$Delta <- (DeltaTab$mexpr - DeltaTab$fexpr)/(DeltaTab$mexpr + DeltaTab$fexpr)

FSTdata <- read.table(paste0(anc,"chr",chrnum,".INFO"),header = TRUE)

annot <- read.table("gencode.v21.annotation.gtf.gz",header = FALSE,
                    stringsAsFactors = FALSE,sep = "\t")

names(annot) <- c("seqname","source","feature","start","end","score","strand",
                  "frame","attribute")
annot <- annot[which(annot$feature == "gene"),]
annot <- annot[which(annot$seqname == paste0("chr",chrnum)),]

annot$Name <- gsub(".*gene_id (.+?);.*","\\1",annot$attribute)
annot <- annot[which(gsub(".*gene_type (.+?);.*","\\1",annot$attribute) == "protein_coding"),]

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

AFXY <- paste0("AF_",tolower(anc),"_XY")
AFXX <- paste0("AF_",tolower(anc),"_XX")
ANXY <- paste0("AN_",tolower(anc),"_XY")
ANXX <- paste0("AN_",tolower(anc),"_XX")
ACXY <- paste0("AC_",tolower(anc),"_XY")
ACXX <- paste0("AC_",tolower(anc),"_XX")

FSTdata[[AFXY]] <- as.numeric(FSTdata[[AFXY]])
FSTdata[[AFXX]] <- as.numeric(FSTdata[[AFXX]])
FSTdata[[ANXY]] <- as.numeric(FSTdata[[ANXY]])
FSTdata[[ANXX]] <- as.numeric(FSTdata[[ANXX]])

FSTdata <- FSTdata[-which(is.na(FSTdata[[AFXY]])),]
FSTdata <- FSTdata[-which(is.na(FSTdata[[AFXX]])),]

AFm <- FSTdata[[AFXY]]
AFf <- FSTdata[[AFXX]]
ANm <- FSTdata[[ANXY]]
ANf <- FSTdata[[ANXX]]

FSTdata$NUMER <- getnumerHud(AFm,AFf,ANm,ANf)
FSTdata$DENOM <- getdenomHud(AFm,AFf)
FSTdata$FST <- FSTdata$NUMER/FSTdata$DENOM

FSTdata <- FSTdata[-which(is.na(FSTdata$FST)),]

FSTdata <- FSTdata[-which((FSTdata[[ACXY]] == 0 & FSTdata[[ACXX]] == 1) | (FSTdata[[ACXY]] == 1 & FSTdata[[ACXX]] == 0)),]

GetGeneRange <- function(genename){
  genestart <- (annot$start[which(annot$Name == genename)] - 1000)
  geneend <- annot$end[which(annot$Name == genename)]
  genechr <- annot$seqname[which(annot$Name == genename)]
  posrangel <- length(which(FSTdata$POS >= genestart & 
                              FSTdata$POS <= geneend &
                              FSTdata$CHROM == genechr))
  if(posrangel == 0){return(0)}
  posrange <- which(FSTdata$POS >= genestart & 
                      FSTdata$POS <= geneend &
                      FSTdata$CHROM == genechr)
  return(posrange)
}
GetGeneAvgFST <- function(posrange){
  avgFST <- sum(FSTdata$NUMER[posrange])/sum(FSTdata$DENOM[posrange])
  return(avgFST)
}
GetGeneMeanFST <- function(posrange){
  meanFST <- mean(FSTdata$FST[posrange])
  return(meanFST)
}
GetGeneMaxFST <- function(posrange){
  maxFST <- max(FSTdata$FST[posrange])
  return(maxFST)
}

chrnumGenes <- annot$Name[which(annot$seqname == paste0("chr",chrnum))]
TwinPeaksDF <- data.frame(Name = DeltaTab$Name[which(DeltaTab$Name %in% chrnumGenes)],
                          Delta = DeltaTab$Delta[which(DeltaTab$Name %in% chrnumGenes)])

TwinPeaksDF$AvgFST <- rep(NA,nrow(TwinPeaksDF))
TwinPeaksDF$MeanFST <- rep(NA,nrow(TwinPeaksDF))
TwinPeaksDF$MaxFST <- rep(NA,nrow(TwinPeaksDF))
TwinPeaksDF$LociCounts <- rep(NA,nrow(TwinPeaksDF))

for(i in 1:nrow(TwinPeaksDF)){
  # print(i)
  # cat("\r",i,"/",nrow(TwinPeaksDF.rand.chr1),sep = "")
  rangei <- GetGeneRange(TwinPeaksDF$Name[i])
  # print(length(rangei))
  if(rangei[1] == 0){next}
  # print("moving on")
  TwinPeaksDF$AvgFST[i] <- GetGeneAvgFST(rangei)
  TwinPeaksDF$MeanFST[i] <- GetGeneMeanFST(rangei)
  TwinPeaksDF$MaxFST[i] <- GetGeneMaxFST(rangei)
  TwinPeaksDF$LociCounts[i] <- length(rangei)
}

if(length(which(is.na(TwinPeaksDF$AvgFST))) > 0){
  TwinPeaksDF <- TwinPeaksDF[-which(is.na(TwinPeaksDF$AvgFST)),]
}


write.table(TwinPeaksDF,file = paste0("./ancs/",anc,"/",samptype,"/TwinPeaksDF.gnomADv3_",anc,".GTExv8.chr",chrnum,".txt"),
            col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t",append = FALSE)

