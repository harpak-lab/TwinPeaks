# NewAnalysis_eQTLRegression_MH.R

library(vroom)

tissue <- "gonads"
pcut <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

# Getting Files

print("Getting Annotation")

annot <- vroom("gencode.v26.annotation.protein_coding.txt")

print("Getting FST Data")

FSTdata <- vroom("./NFEchrALL.cisFilt.INFO")

print("Getting Testis eQTL Data")

eQTLlist_testis <- vroom(paste0("./eQTLs/testis.filtered.MHcorrection.NEW.txt"))
eQTLlist_testis$BONAll <- p.adjust(eQTLlist_testis$pvalue,method = "bonferroni")
eQTLlist_testis <- eQTLlist_testis[which(eQTLlist_testis$BONAll < pcut),]

print("Getting Ovary eQTL Data")

eQTLlist_ovary <- vroom(paste0("./eQTLs/ovary.filtered.MHcorrection.NEW.txt"))
eQTLlist_ovary$BONAll <- p.adjust(eQTLlist_ovary$pvalue,method = "bonferroni")
eQTLlist_ovary <- eQTLlist_ovary[which(eQTLlist_ovary$BONAll < pcut),]

print("Filtering FST Data")

FSTdata <- FSTdata[which(FSTdata$Loci %in% c(unique(eQTLlist_testis$Loci),unique(eQTLlist_ovary$Loci))),]
eQTLlist_testis <- eQTLlist_testis[which(eQTLlist_testis$Loci %in% FSTdata$Loci),]
eQTLlist_ovary <- eQTLlist_ovary[which(eQTLlist_ovary$Loci %in% FSTdata$Loci),]

eQTLlist <- rbind(eQTLlist_testis,eQTLlist_ovary)

FSTdata$FSTSE <- rep(NA,nrow(FSTdata))
for(i in 1:nrow(FSTdata)){
  set.seed(i)
  print(paste0(i,"/",nrow(FSTdata)))
  pm_boot <- rbinom(10000,FSTdata$AN_nfe_XY[i],FSTdata$AF_nfe_XY[i])
  pf_boot <- rbinom(10000,FSTdata$AN_nfe_XX[i],FSTdata$AF_nfe_XX[i])
  nm_boot <- rep(FSTdata$AN_nfe_XY[i],10000)
  nf_boot <- rep(FSTdata$AN_nfe_XX[i],10000)
  FSTboot <- (((pm_boot-pf_boot)^2) - ((pm_boot*(1-pm_boot))/(nm_boot-1)) - ((pf_boot*(1-pf_boot))/(nf_boot-1)))/((pm_boot*(1-pf_boot)) + (pf_boot*(1-pm_boot)))
  FSTdata$FSTSE[i] <- sd(FSTboot)
}

print("Getting Tissue Expression Data")

testisTPM <- vroom(paste0("GTExTPM/gene_tpm_2017-06-05_v8_testis.gct.gz"),
                delim = "\t",skip = 2)
ovaryTPM <- vroom(paste0("GTExTPM/gene_tpm_2017-06-05_v8_ovary.gct.gz"),
                   delim = "\t",skip = 2)

testisTPM$NameAbr <- gsub("\\..*","",testisTPM$Name)
ovaryTPM$NameAbr <- gsub("\\..*","",ovaryTPM$Name)

testisTPM <- testisTPM[which(testisTPM$NameAbr %in% ovaryTPM$NameAbr),]
ovaryTPM <- ovaryTPM[which(ovaryTPM$NameAbr %in% testisTPM$NameAbr),]

tisTPM <- cbind(testisTPM[,c(1:3,ncol(testisTPM),4:(ncol(testisTPM)-1))],ovaryTPM[,c(4:(ncol(ovaryTPM)-1))])

print("Filtering Expr Data")

# tisTPM$NameAbr <- gsub("\\..*","",tisTPM$Name)
tisTPM <- tisTPM[which(tisTPM$NameAbr %in% eQTLlist$gene_id),]

sampleAtt <- vroom("GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")
Males <- sampleAtt$SUBJID[which(sampleAtt$SEX == 1)]
Females <- sampleAtt$SUBJID[which(sampleAtt$SEX == 2)]

tisTPMMales <- grep(paste(Males,collapse = "|"),names(tisTPM))
tisTPMFemales <- grep(paste(Females,collapse = "|"),names(tisTPM))
tisTPMAll <- c(tisTPMMales,tisTPMFemales)

print("Getting Delta Values")

tisTPM$MAvg <- rowMeans(tisTPM[1:nrow(tisTPM),tisTPMMales])
tisTPM$FAvg <- rowMeans(tisTPM[1:nrow(tisTPM),tisTPMFemales])
tisTPM$Delta <- (tisTPM$MAvg - tisTPM$FAvg)/(tisTPM$MAvg + tisTPM$FAvg)

DeltaSEMat <- matrix(NA,ncol = 10000,nrow = nrow(tisTPM))
for(i in 1:10000){
  set.seed(i)
  print(paste0(i,"/10000"))
  bootMales <- sample(tisTPMMales,replace = TRUE)
  bootFemales <- sample(tisTPMFemales,replace = TRUE)
  m_mean_i <- rowMeans(tisTPM[,bootMales])
  f_mean_i <- rowMeans(tisTPM[,bootFemales])
  DeltaSEMat[,i] <- (m_mean_i - f_mean_i)/(m_mean_i + f_mean_i)
}

RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

DeltaSE2 <- RowVar(DeltaSEMat)

print("Getting Regression 1")

RegressionDF <- data.frame(Gene = tisTPM$NameAbr,
                           Delta = tisTPM$Delta,
                           DeltaSE = DeltaSE2)

RegressionDF$FST <- rep(NA,nrow(RegressionDF))
RegressionDF$PTOT <- rep(NA,nrow(RegressionDF))
RegressionDF$FSTSE <- rep(NA,nrow(RegressionDF))

OutDF <- data.frame(AVal = rep(NA,1),
                    ASE = rep(NA,1),
                    Int = rep(NA,1))

sma_weighted <- function(Delta, p, y, sx, sy) {
  
  #remove NAs
  NArows <- union(union(union(which(is.na(Delta)),which(is.na(p))),which(is.na(y))),union(which(is.na(sx)),which(is.na(sy))))
  
  if(length(NArows) > 0){
    Delta <- Delta[-NArows]
    p <- p[-NArows]
    y <- y[-NArows]
    sx <- sx[-NArows]
    sy <- sy[-NArows]
  }
  
  x <- 4*p*(1-p)*(Delta^2)
  
  #weights
  wx <- sx^2
  wy <- sy^2
  hm <- 2/((1/wx) + (1/wy))
  v <- 1/hm
  w <- v/sum(v)
  
  #weighted means  
  x_bar <- sum(w * x)
  y_bar <- sum(w * y)
  
  #weighted covariance and variances
  cov_xy <- sum(w * (x - x_bar) * (y - y_bar))
  var_x <- sum(w * (x - x_bar)^2)
  var_y <- sum(w * (y - y_bar)^2)
  
  #slope and intercept
  b1 <- sign(cov_xy)*(var_y / var_x)^0.5
  b0 <- y_bar - b1 * x_bar
  
  #SE
  n <- length(x)
  r <- cov_xy/sqrt(var_x*var_y)
  se<- abs(b1) * sqrt((1-r^2)/n)
  
  return(c(b0, b1,se))
}

for(i in 1:nrow(RegressionDF)){
  print(paste0(i,"/",nrow(RegressionDF)))
  # print(i)
  # cat(i,"/",nrow(RegDF),"\r")
  gene_i <- RegressionDF$Gene[i]
  # CisLoci <- eQTLlist$Loci[which(eQTLlist$chromosome == annot$seqname[which(annot$NameAbr == gene_i)] &
  #                                  eQTLlist$position < (annot$midpoint[which(annot$NameAbr == gene_i)] + 1000000) &
  #                                  eQTLlist$position < (annot$midpoint[which(annot$NameAbr == gene_i)] - 1000000))]
  BestLocus <- eQTLlist$Loci[which(eQTLlist$pvalue == min(eQTLlist$pvalue[which(eQTLlist$gene_id == gene_i)]))]
  if(length(BestLocus > 1)){
    PossibleLoci <- data.frame(Loci = FSTdata$Loci[which(FSTdata$Loci %in% BestLocus)],
                               PTOT = FSTdata$PTOT[which(FSTdata$Loci %in% BestLocus)])
    BestLocus <- PossibleLoci$Loci[sample(which(PossibleLoci$PTOT == min(PossibleLoci$PTOT)),1)]
  }
  RegressionDF$FST[i] <- FSTdata$FST[which(FSTdata$Loci == BestLocus)]
  RegressionDF$PTOT[i] <- FSTdata$PTOT[which(FSTdata$Loci == BestLocus)]
  RegressionDF$FSTSE[i] <- FSTdata$FSTSE[which(FSTdata$Loci == BestLocus)]
}

if(length(which(is.na(RegressionDF$FST))) > 0){RegressionDF <- RegressionDF[-which(is.na(RegressionDF$FST)),]}

outvec <- sma_weighted(RegressionDF$Delta,RegressionDF$PTOT,RegressionDF$FST,
                       RegressionDF$DeltaSE,RegressionDF$FSTSE)

OutDF$AVal[1] <- outvec[2]
OutDF$ASE[1] <- outvec[3]
OutDF$Int[1] <- outvec[1]
# OutDF$IntSE[1] <- summary(AModel)$coefficients[1,2]

write.table(OutDF,file = paste0("./AVals_New_MH_Rev4/AVals_NFE.PVal",pcut,".",tissue,".MHNew.Revision.txt"),
            col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t",append = FALSE)

print("Iterating Regression")

set.seed(1)
tisTPMAll <- sample(tisTPMAll)
tisTPMMales <- tisTPMAll[1:length(tisTPMMales)]
tisTPMFemales <- tisTPMAll[(length(tisTPMMales) + 1):length(tisTPMAll)]

tisTPM$MAvg <- rowMeans(tisTPM[1:nrow(tisTPM),tisTPMMales])
tisTPM$FAvg <- rowMeans(tisTPM[1:nrow(tisTPM),tisTPMFemales])
tisTPM$Delta <- (tisTPM$MAvg - tisTPM$FAvg)/(tisTPM$MAvg + tisTPM$FAvg)

RegressionDF$Delta <- tisTPM$Delta

for(i in 1:10000){
  set.seed(i)
  print(paste0(i,"/10000"))
  bootMales <- sample(tisTPMMales,replace = TRUE)
  bootFemales <- sample(tisTPMFemales,replace = TRUE)
  m_mean_i <- rowMeans(tisTPM[,bootMales])
  f_mean_i <- rowMeans(tisTPM[,bootFemales])
  DeltaSEMat[,i] <- (m_mean_i - f_mean_i)/(m_mean_i + f_mean_i)
}

DeltaSE2 <- RowVar(DeltaSEMat)
RegressionDF$DeltaSE <- DeltaSE2

outvec <- sma_weighted(RegressionDF$Delta,RegressionDF$PTOT,RegressionDF$FST,
                       RegressionDF$DeltaSE,RegressionDF$FSTSE)

OutDF$AVal[1] <- outvec[2]
OutDF$ASE[1] <- outvec[3]
OutDF$Int[1] <- outvec[1]

write.table(OutDF,file = paste0("./AVals_New_MH_Rev4/AVals_NFE.PVal",pcut,".",tissue,".MHNew.Revision.txt"),
            col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t",append = TRUE)

for(i in 2:1000){
  set.seed(i)
  print(paste0(i,"/1000"))
  tisTPMAll <- sample(tisTPMAll)
  tisTPMMales <- tisTPMAll[1:length(tisTPMMales)]
  tisTPMFemales <- tisTPMAll[(length(tisTPMMales) + 1):length(tisTPMAll)]
  
  tisTPM$MAvg <- rowMeans(tisTPM[1:nrow(tisTPM),tisTPMMales])
  tisTPM$FAvg <- rowMeans(tisTPM[1:nrow(tisTPM),tisTPMFemales])
  tisTPM$Delta <- (tisTPM$MAvg - tisTPM$FAvg)/(tisTPM$MAvg + tisTPM$FAvg)

  RegressionDF$Delta <- tisTPM$Delta

  # AModel <- lm(FST ~ XVar,data = RegressionDF,weights = 1/Var)
  
  outvec <- sma_weighted(RegressionDF$Delta,RegressionDF$PTOT,RegressionDF$FST,
                         RegressionDF$DeltaSE,RegressionDF$FSTSE)
  
  OutDF$AVal[1] <- outvec[2]
  OutDF$ASE[1] <- outvec[3]
  OutDF$Int[1] <- outvec[1]
  
  write.table(OutDF,file = paste0("./AVals_New_MH_Rev4/AVals_NFE.PVal",pcut,".",tissue,".MHNew.Revision.txt"),
              col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t",append = TRUE)
}

# To run
# Rscript >AVals_New_MH_Rev4/gonads_AVal.05.Bonf.Rev.out 2>AVals_New_MH_Rev4/gonads_AVal.05.Bonf.Rev.err NewAnalysis_eQTLRegression_MH_Gonads_Revision4.R "0.05" &
