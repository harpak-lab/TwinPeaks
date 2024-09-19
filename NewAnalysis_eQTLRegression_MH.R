# NewAnalysis_eQTLRegression_MH.R

library(vroom)

tissue <- commandArgs(trailingOnly = TRUE)[1]
pcut <- as.numeric(commandArgs(trailingOnly = TRUE)[2])

# Getting Files

print("Getting Annotation")

annot <- vroom("gencode.v26.annotation.protein_coding.txt")

print("Getting FST Data")

FSTdata <- vroom("./NFEchrALL.cisFilt.INFO")

print("Getting eQTL Data")

eQTLlist <- vroom(paste0("./eQTLs/",tissue,".filtered.MHcorrection.txt"))
eQTLlist$BONAll <- p.adjust(eQTLlist$pvalue,method = "bonferroni")
eQTLlist <- eQTLlist[which(eQTLlist$BONAll < pcut),]

print("Filtering FST Data")

FSTdata <- FSTdata[which(FSTdata$Loci %in% unique(eQTLlist$Loci)),]
eQTLlist <- eQTLlist[which(eQTLlist$Loci %in% FSTdata$Loci),]

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
tisTPM$Var <- apply(tisTPM[1:nrow(tisTPM),c(tisTPMMales,tisTPMFemales)],1,var)

print("Getting Regression 1")

RegressionDF <- data.frame(Gene = tisTPM$NameAbr,
                           Delta = tisTPM$Delta,
                           Var = tisTPM$Var)

RegressionDF$FST <- rep(NA,nrow(RegressionDF))
RegressionDF$PTOT <- rep(NA,nrow(RegressionDF))

OutDF <- data.frame(AVal = rep(NA,1),
                    ASE = rep(NA,1),
                    Int = rep(NA,1),
                    IntSE = rep(NA,1))

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
}

RegressionDF$XVar <- 4*RegressionDF$PTOT*(1-RegressionDF$PTOT)*(RegressionDF$Delta^2)

AModel <- lm(FST ~ XVar,data = RegressionDF,weights = 1/Var)

OutDF$AVal[1] <- summary(AModel)$coefficients[2,1]
OutDF$ASE[1] <- summary(AModel)$coefficients[2,2]
OutDF$Int[1] <- summary(AModel)$coefficients[1,1]
OutDF$IntSE[1] <- summary(AModel)$coefficients[1,2]

write.table(OutDF,file = paste0("./AVals_New_MH/AVals_NFE.PVal",pcut,".",tissue,".MHNew.txt"),
            col.names = TRUE,row.names = FALSE,quote = FALSE,sep = "\t",append = FALSE)

print("Iterating Regression")

for(i in 2:1000){
  print(paste0(i,"/1000"))
  tisTPMAll <- sample(tisTPMAll)
  tisTPMMales <- tisTPMAll[1:length(tisTPMMales)]
  tisTPMFemales <- tisTPMAll[(length(tisTPMMales) + 1):length(tisTPMAll)]
  
  tisTPM$MAvg <- rowMeans(tisTPM[1:nrow(tisTPM),tisTPMMales])
  tisTPM$FAvg <- rowMeans(tisTPM[1:nrow(tisTPM),tisTPMFemales])
  tisTPM$Delta <- (tisTPM$MAvg - tisTPM$FAvg)/(tisTPM$MAvg + tisTPM$FAvg)
  tisTPM$Var <- apply(tisTPM[1:nrow(tisTPM),c(tisTPMMales,tisTPMFemales)],1,var)
  
  RegressionDF$Delta <- tisTPM$Delta
  RegressionDF$Var <- tisTPM$Var
  RegressionDF$XVar <- 4*RegressionDF$PTOT*(1-RegressionDF$PTOT)*(RegressionDF$Delta^2)
  
  AModel <- lm(FST ~ XVar,data = RegressionDF,weights = 1/Var)
  
  OutDF$AVal[1] <- summary(AModel)$coefficients[2,1]
  OutDF$ASE[1] <- summary(AModel)$coefficients[2,2]
  OutDF$Int[1] <- summary(AModel)$coefficients[1,1]
  OutDF$IntSE[1] <- summary(AModel)$coefficients[1,2]
  
  write.table(OutDF,file = paste0("./AVals_New_MH/AVals_NFE.PVal",pcut,".",tissue,".MHNew.txt"),
              col.names = FALSE,row.names = FALSE,quote = FALSE,sep = "\t",append = TRUE)
}



# To run
# Rscript >AVals_New_MH/adipose_subcutaneous_AVal.05.Bonf.out 2>AVals_New_MH/adipose_subcutaneous_AVal.05.Bonf.err NewAnalysis_eQTLRegression_MH.R adipose_subcutaneous "0.05" &
# for tis in artery_aorta brain_cerebellar_hemisphere brain_cortex brain_hippocampus brain_spinal_cord_cervical_c-1;do echo ${tis};(Rscript >AVals_New_MH/${tis}_AVal.05.Bonf.out 2>AVals_New_MH/${tis}_AVal.05.Bonf.err NewAnalysis_eQTLRegression_MH.R ${tis} "0.05" &);done
# for tis in cells_ebv-transformed_lymphocytes esophagus_mucosa kidney_cortex muscle_skeletal nerve_tibial;do echo ${tis};(Rscript >AVals_New_MH/${tis}_AVal.05.Bonf.out 2>AVals_New_MH/${tis}_AVal.05.Bonf.err NewAnalysis_eQTLRegression_MH.R ${tis} "0.05" &);done
# for tis in skin_sun_exposed_lower_leg stomach whole_blood;do echo ${tis};(Rscript >AVals_New_MH/${tis}_AVal.05.Bonf.out 2>AVals_New_MH/${tis}_AVal.05.Bonf.err NewAnalysis_eQTLRegression_MH.R ${tis} "0.05" &);done

