# NewAnalysis_eQTLRegression_MH.R

library(vroom)

tissue <- commandArgs(trailingOnly = TRUE)[1]
pcut <- as.numeric(commandArgs(trailingOnly = TRUE)[2])

# tissue <- "adipose_subcutaneous"
# pcut <- 0.05

# Getting Files

print("Getting Annotation")

annot <- vroom("gencode.v26.annotation.protein_coding.txt")

print("Getting FST Data")

FSTdata <- vroom("./NFEchrALL.cisFilt.INFO")

print("Getting eQTL Data")

eQTLlist <- vroom(paste0("./eQTLs/",tissue,".filtered.MHcorrection.NEW.txt"))
eQTLlist$BONAll <- p.adjust(eQTLlist$pvalue,method = "bonferroni")
eQTLlist <- eQTLlist[which(eQTLlist$BONAll < pcut),]

print("Filtering FST Data")

FSTdata <- FSTdata[which(FSTdata$Loci %in% unique(eQTLlist$Loci)),]
eQTLlist <- eQTLlist[which(eQTLlist$Loci %in% FSTdata$Loci),]

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

tisTPM <- vroom(paste0("GTExTPM/gene_tpm_2017-06-05_v8_",tissue,".gct.gz"),
                delim = "\t",skip = 2)

print("Filtering Expr Data")

tisTPM$NameAbr <- gsub("\\..*","",tisTPM$Name)
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
# Rscript >AVals_New_MH_Rev4/adipose_subcutaneous_AVal.05.Bonf.Rev.out 2>AVals_New_MH_Rev4/adipose_subcutaneous_AVal.05.Bonf.Rev.err NewAnalysis_eQTLRegression_MH_Revision4.R adipose_subcutaneous "0.05" &
# for tis in esophagus_muscularis heart_atrial_appendage heart_left_ventricle colon_transverse esophagus_gastroesophageal_junction;do echo ${tis};(Rscript >AVals_New_MH_Rev4/${tis}_AVal.05.Bonf.Rev.out 2>AVals_New_MH_Rev4/${tis}_AVal.05.Bonf.Rev.err NewAnalysis_eQTLRegression_MH_Revision4.R ${tis} "0.05" &);done
# for tis in pituitary skin_not_sun_exposed_suprapubic small_intestine_terminal_ileum spleen thyroid;do echo ${tis};(Rscript >AVals_New_MH_Rev4/${tis}_AVal.05.Bonf.Rev.out 2>AVals_New_MH_Rev4/${tis}_AVal.05.Bonf.Rev.err NewAnalysis_eQTLRegression_MH_Revision4.R ${tis} "0.05" &);done
# for tis in adipose_visceral_omentum adrenal_gland artery_coronary artery_tibial brain_amygdala brain_anterior_cingulate_cortex_ba24;do echo ${tis};(Rscript >AVals_New_MH_Rev4/${tis}_AVal.05.Bonf.Rev.out 2>AVals_New_MH_Rev4/${tis}_AVal.05.Bonf.Rev.err NewAnalysis_eQTLRegression_MH_Revision4.R ${tis} "0.05" &);done
# for tis in brain_caudate_basal_ganglia brain_cerebellum brain_frontal_cortex_ba9 brain_hypothalamus;do echo ${tis};(Rscript >AVals_New_MH_Rev4/${tis}_AVal.05.Bonf.Rev.out 2>AVals_New_MH_Rev4/${tis}_AVal.05.Bonf.Rev.err NewAnalysis_eQTLRegression_MH_Revision4.R ${tis} "0.05" &);done
# for tis in brain_nucleus_accumbens_basal_ganglia brain_putamen_basal_ganglia brain_substantia_nigra breast_mammary_tissue cells_cultured_fibroblasts colon_sigmoid;do echo ${tis};(Rscript >AVals_New_MH_Rev4/${tis}_AVal.05.Bonf.Rev.out 2>AVals_New_MH_Rev4/${tis}_AVal.05.Bonf.Rev.err NewAnalysis_eQTLRegression_MH_Revision4.R ${tis} "0.05" &);done
# for tis in artery_aorta brain_cerebellar_hemisphere brain_cortex brain_hippocampus;do echo ${tis};(Rscript >AVals_New_MH_Rev4/${tis}_AVal.05.Bonf.Rev.out 2>AVals_New_MH_Rev4/${tis}_AVal.05.Bonf.Rev.err NewAnalysis_eQTLRegression_MH_Revision4.R ${tis} "0.05" &);done
# for tis in brain_spinal_cord_cervical_c-1 cells_ebv-transformed_lymphocytes esophagus_mucosa kidney_cortex muscle_skeletal;do echo ${tis};(Rscript >AVals_New_MH_Rev4/${tis}_AVal.05.Bonf.Rev.out 2>AVals_New_MH_Rev4/${tis}_AVal.05.Bonf.Rev.err NewAnalysis_eQTLRegression_MH_Revision4.R ${tis} "0.05" &);done
# for tis in nerve_tibial skin_sun_exposed_lower_leg stomach whole_blood;do echo ${tis};(Rscript >AVals_New_MH_Rev4/${tis}_AVal.05.Bonf.Rev.out 2>AVals_New_MH_Rev4/${tis}_AVal.05.Bonf.Rev.err NewAnalysis_eQTLRegression_MH_Revision4.R ${tis} "0.05" &);done
