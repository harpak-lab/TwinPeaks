# Twin Peaks Figs

library(ggplot2)
library(ggpubr)
library(vroom)
library(ggrepel)
library(grid)
library(ggstance)
# library(pBrackets)

# Cartoon Fig

xvals <- seq(-1,1,by = 0.01)
yvals <- (xvals - 1)*(xvals + 1)*(xvals^2)*(-1/5)
cartoondf <- data.frame(Delta = xvals,
                        FST = yvals)
cartoondfLow <- cartoondf[which(cartoondf$Delta >= -0.4 & cartoondf$Delta <= 0.4),]
cartoonfig <- ggplot() +
  geom_hline(yintercept = 0,linetype = "dashed") + geom_vline(xintercept = 0,linetype = "dashed") +
  geom_line(data = cartoondf,mapping = aes(x = Delta,y = FST),linetype = "dashed",linewidth = 1.25,color = "grey") +
  geom_point(aes(x = -1,y = 0),color = "grey",size = 6) +
  geom_point(aes(x = 1,y = 0),color = "grey",size = 6) +
  geom_line(data = cartoondfLow,mapping = aes(x = Delta,y = FST),linewidth = 1.25,color = "blue") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",legend.title = element_blank(),
        text = element_text(size = 14)) +
  # ggtitle("Theoretical expectation for Twin Peaks pattern is unsubstantiated") +
  xlab("Sex difference in gene expression, \u0394") +
  ylab(expression("Between-Sex "~F[ST])) +
  coord_cartesian(xlim = c(-1,1),ylim = c(-0.01,0.05)) +
  scale_x_continuous(breaks = c(-1,0,1)) + scale_y_continuous(breaks = c(0)) +
  geom_text_repel(aes(x = -1,y = 0),color = "black",label = "Unsupported\nIntuition",
                  nudge_x = 0.3,nudge_y = 0.005,size = 4) +
  geom_text_repel(aes(x = 1,y = 0),color = "black",label = "Unsupported\nIntuition",
                  nudge_x = -0.3,nudge_y = 0.005,size = 4) +
  # geom_text_repel(aes(x = -1,y = 0),color = "black",label = "Interpolation",
  #                 nudge_x = 0.3,nudge_y = 0.025) +
  geom_text_repel(aes(x = -0.9,y = 0.03078),color = "black",label = "Interpolation",
                  nudge_x = 0.2,nudge_y = -0.005,size = 4) +
  # geom_text_repel(aes(x = 1,y = 0),color = "black",label = "Interpolation",
  #                 nudge_x = -0.3,nudge_y = 0.025) +
  geom_text_repel(aes(x = 0.9,y = 0.03078),color = "black",label = "Interpolation",
                  nudge_x = -0.2,nudge_y = -0.005,size = 4) +
  geom_label_repel(aes(x = -0.3,y = 0.01638),color = "blue",label = "Theoretically Supported Model",
                  nudge_x = 0.3,nudge_y = 0.02362,size = 4) +
  geom_label_repel(aes(x = 0.3,y = 0.01638),color = "blue",label = "Theoretically Supported Model",
                  nudge_x = -0.3,nudge_y = 0.02362,size = 4)

# png("./RevFigs/CartoonFig.png",width = 750,height = 450)
# print(cartoonfig)
# dev.off()
# 
# ggsave("./RevFigs/CartoonFig.svg",plot = cartoonfig,
#        width = 3000,height = 1800,units = "px")

# ggsave("./RevFigs/Fig1.png",plot = cartoonfig,
       # width = 3000,height = 1800,units = "px")

# ggsave("./RevFigs/Fig1.pdf",plot = cartoonfig,
#        width = 3000,height = 1800,units = "px",device = cairo_pdf)

ggsave("./RevFigs/Fig1_FINAL.pdf",plot = cartoonfig,
       width = 6.5,height = 3.9,units = "in",dpi = 300,device = cairo_pdf)


# svg("./RevFigs/CartoonFig.svg",
#     width = 750,height = 450)
# print(cartoonfig)
# dev.off()

# Getting dataframes
randTP_GTEx <- read.table("AllAncTP.RAND_GTEx.avgHUD.All.txt",
                          header = TRUE,stringsAsFactors = FALSE)
norm_counts <- read.table("norm_counts_filt.Ensembl77.GTEx_v3.txt",
                          header = TRUE,stringsAsFactors = FALSE)

# Getting coefficients of each iteration
mOG <- lm(randTP_GTEx$FST ~ poly(randTP_GTEx$deltaOG,4,raw = TRUE))
coef_df <- data.frame(X4 = rep(NA,500),
                      X3 = rep(NA,500),
                      X2 = rep(NA,500),
                      X1 = rep(NA,500),
                      Int = rep(NA,500))

for(i in 1:500){
  if(i %% 50 == 0){print(i)}
  # print(i)
  colnum <- which(names(randTP_GTEx) == paste0("delta",i))
  y <- randTP_GTEx$FST
  x <- randTP_GTEx[[colnum]]
  if(length(which(is.na(x))) > 0){
    y <- y[-which(is.na(x))]
    x <- x[-which(is.na(x))]
  }
  # print(length(x))
  # print(length(y))
  modeli <- lm(y ~ poly(x,4,raw = TRUE))
  coef_df$X4[i] <- coef(modeli)[5]
  coef_df$X3[i] <- coef(modeli)[4]
  coef_df$X2[i] <- coef(modeli)[3]
  coef_df$X1[i] <- coef(modeli)[2]
  coef_df$Int[i] <- coef(modeli)[1]
}

# Supp 1. Coefficients of Permuted Iterations

# pdf("TwinPeaksFigs12.pdf")
# ggplot(coef_df) + geom_histogram(aes(X4)) + geom_vline(xintercept = coef(mOG)[5])
# ggplot(coef_df) + geom_histogram(aes(X3)) + geom_vline(xintercept = coef(mOG)[4])
# ggplot(coef_df) + geom_histogram(aes(X2)) + geom_vline(xintercept = coef(mOG)[3])
# ggplot(coef_df) + geom_histogram(aes(X1)) + geom_vline(xintercept = coef(mOG)[2])
# ggplot(coef_df) + geom_histogram(aes(Int),binwidth = 1e-15) + geom_vline(xintercept = coef(mOG)[1])

# length(which(coef_df$X4 < coef(mOG)[5]))/500
# length(which(coef_df$X3 < coef(mOG)[4]))/500
# length(which(coef_df$X2 < coef(mOG)[3]))/500
# length(which(coef_df$X1 < coef(mOG)[2]))/500
# length(which(coef_df$Int < coef(mOG)[1]))/500

# Checking the C&K criteria for Twin Peaks

CKTest <- data.frame(ANOVA = rep(NA,500),
                     X4Coeff = coef_df$X4,
                     RealRoots = rep(NA,500),
                     PassANOVA = rep(FALSE,500),
                     PassCoeff = rep(FALSE,500),
                     PassRoots = rep(FALSE,500),
                     PassTest = rep(FALSE,500))

# ANOVA Test

ANOVATest <- function(i){
  colnum <- paste0("delta",i)
  y <- randTP_GTEx$FST
  x <- randTP_GTEx[[colnum]]
  if(length(which(is.na(x))) > 0){
    y <- y[-which(is.na(x))]
    x <- x[-which(is.na(x))]
  }
  
  m4 <- lm(y ~ poly(x,4,raw = TRUE))
  m3 <- lm(y ~ poly(x,3,raw = TRUE))
  
  pval <- anova(m4,m3)$"Pr(>F)"[2]
  return(pval)
}
for(i in 1:500){
  if(i %% 50 == 0){print(i)}
  CKTest$ANOVA[i] <- ANOVATest(i)
}

CKTest$PassANOVA[which(CKTest$ANOVA < 0.05)] <- TRUE

# X4 Coefficients Test

CKTest$X4Coeff <- coef_df$X4
CKTest$PassCoeff[which(CKTest$X4Coeff < 0)] <- TRUE

# Derivative Roots Test

derivative_coef_df <- data.frame(X3 = 4*coef_df$X4,
                                 X2 = 3*coef_df$X3,
                                 X1 = 2*coef_df$X2,
                                 Int = coef_df$X1)

NumberRealRoots <- function(i){
  roots <- polyroot(z = c(derivative_coef_df[i,4],derivative_coef_df[i,3],
                          derivative_coef_df[i,2],derivative_coef_df[i,1]))
  #Arbitrarily set the cut-off of the imaginary term at 1e-12
  ImCoeff <- Im(roots)
  RealRoots <- roots[which(abs(ImCoeff) < 1e-12)]
  return(length(RealRoots))
}
for(i in 1:500){
  if(i %% 50 == 0){print(i)}
  CKTest$RealRoots[i] <- NumberRealRoots(i)
}

CKTest$PassRoots[which(CKTest$RealRoots == 3)] <- TRUE

CKTest$PassTest[which(CKTest$PassANOVA & CKTest$PassCoeff & CKTest$PassRoots)] <- TRUE

# Fig 1a. FST vs Expr Variance

getTrendPointsFST <- function(sex,n){
  dAll <- data.frame(x = exprVar[[paste0("Var",sex)]],y = exprVar$FST)
  dAll <- dAll[order(dAll$x),]
  binstart <- round(seq(1,17523,length.out = (n+1)))[1:n]
  binend <- round(seq(1,17523,length.out = (n+1)))[2:(n+1)] - 1
  VarTrends <- data.frame(x = rep(NA,n),
                          y = rep(NA,n))
  for(i in 1:n){
    starti <- binstart[i]
    endi <- binend[i]
    VarTrends$x[i] <- mean(dAll$x[starti:endi])
    VarTrends$y[i] <- mean(dAll$y[starti:endi])
  }
  if(length(which(is.infinite(VarTrends$x))) > 0){
    VarTrends <- VarTrends[-which(is.infinite(VarTrends$x)),]
  }
  return(VarTrends)
}

length(which(randTP_GTEx$name %in% norm_counts$GeneName))
norm_counts <- norm_counts[which(norm_counts$GeneName %in% randTP_GTEx$name),]
length(which(norm_counts$GeneName == randTP_GTEx$name))

exprVar <- data.frame(Name = randTP_GTEx$name,
                      FST = randTP_GTEx$FST,
                      Delta = randTP_GTEx$deltaOG)
exprAll <- 3:22
exprFems <- 3:8
exprMales <- 9:22
exprVar$VarAll <- apply(norm_counts[,exprAll],1,var)
exprVar$VarF <- apply(norm_counts[,exprFems],1,var)
exprVar$VarM <- apply(norm_counts[,exprMales],1,var)
exprVar$VarAvg <- (exprVar$VarF + exprVar$VarM)/2
exprVar$VarAll <- log(exprVar$VarAll,base = 10)
exprVar$VarF <- log(exprVar$VarF,base = 10)
exprVar$VarM <- log(exprVar$VarM,base = 10)
exprVar$VarAvg <- log(exprVar$VarAvg,base = 10)

AllTrendsFST <- getTrendPointsFST("Avg",100)
exprvarplot <- ggplot() +
  # geom_point(data = exprVar,aes(x = VarM,y = FST)) +
  geom_point(data = AllTrendsFST,aes(x = x,y = y),color = "blue",size = 1.5) +
  # geom_smooth(data=exprVar,aes(x = VarAll,y = FST),method = "lm",formula = y ~ x,color = "blue",se = FALSE)+
  geom_smooth(data = exprVar,aes(x = VarAvg,y = FST),color = "black",method = "loess",se = FALSE,linewidth = 1) +
  xlab("Sex-averaged variance in gene expression") +
  ylab(expression(atop("Between-Sex ",F[ST]~"(x"*10^{-4}*")"))) +
  ggtitle(expression(F[ST]~" tends to increase with variance in gene expression")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",legend.title = element_blank(),
        plot.title = element_text(margin=margin(t=25,b=-25),
                                  hjust = 0.5),
        text = element_text(size = 9.5),
        plot.margin = unit(c(0, 0.15, 0.1, 0.15),"inches")) +
  # geom_text_repel(aes(x = 4,y = -0.00002),color = "black",label = "LOESS Fit",
  #                 nudge_x = -1,nudge_y = -0.00019) +
  coord_cartesian(xlim = c(-5,5),ylim = c(-0.0003,0.0002)) +
  scale_x_continuous(labels = parse(text = paste0("10^",-5:5)),
                     breaks = -5:5) +
  scale_y_continuous(labels = c(-3,-2,-1,0,1),
                     breaks = c(-0.0003,-0.0002,-0.0001,0,0.0001)) +
  annotation_logticks(sides = "b")

# print(exprvarplot)
# Fig 1b. Binned point distribution

deltabin_GTEx <- function(num,nbins){
  colnum <- which(names(randTP_GTEx) == paste0("delta",num))
  bins <- seq(-1,1,length.out = (nbins + 1))
  Bins.DF <- data.frame(delta = rep(NA,nbins),
                        meanFST = rep(NA,nbins),
                        sdFST = rep(NA,nbins))
  for(i in 1:nbins){
    binstart <- bins[i]
    binend <- bins[i+1]
    Bins.DF$delta[i] <- mean(c(binstart,binend))
    bingenes <- which(randTP_GTEx[[colnum]] >= binstart & randTP_GTEx[[colnum]] <= binend)
    Bins.DF$meanFST[i] <- mean(randTP_GTEx$FST[bingenes])
    Bins.DF$sdFST[i] <- sd(randTP_GTEx$FST[bingenes])
  }
  return(Bins.DF)
}

GTEx_DeltaOGBins40 <- deltabin_GTEx("OG",40)
# GTEx_DeltaOGBins40$delta <- factor(GTEx_DeltaOGBins40$delta)
binplot <- ggplot()
for(i in which(!CKTest$PassTest)){
  # print(i)
  GTEx_DeltaiBins40 <- deltabin_GTEx(i,40)
  # GTEx_DeltaiBins40$delta <- factor(GTEx_DeltaiBins40$delta)
  binplot <- binplot +
    geom_point(data = GTEx_DeltaiBins40,aes(x = delta,y = meanFST,color = "fail"),size = 1)
}
for(i in which(CKTest$PassTest)){
  # print(i)
  GTEx_DeltaiBins40 <- deltabin_GTEx(i,40)
  # GTEx_DeltaiBins40$delta <- factor(GTEx_DeltaiBins40$delta)
  binplot <- binplot +
    geom_point(data = GTEx_DeltaiBins40,aes(x = delta,y = meanFST,color = "pass"),size = 1)
}
binplot <- binplot + geom_point(data = GTEx_DeltaOGBins40,aes(x = delta,y = meanFST,color = "real"),
                                          size = 2) +
  xlab(expression("Sex difference in gene expression, \u0394")) +
  ylab(expression(atop("Between-Sex ",F[ST]~"(x"*10^{-3}*")"))) +
  ggtitle("Effect of sex differential gene expression is consistent with empirical null") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",legend.title = element_blank(),
        legend.background = element_blank(),
        plot.title = element_text(margin=margin(t=25,b=-25),hjust = 0.5),
        text = element_text(size = 9.5),
        plot.margin = unit(c(0, 0.15, 0, 0.15),"inches")) +
  scale_color_manual(values = c(fail = "orange",pass = "red",real = "blue"),
                     labels = c(fail = "Non-Twin Peaks Iterations",
                                pass = "Twin Peaks Iterations",
                                real = "Real Data")) +
  guides(color = guide_legend(nrow = 1)) +
  scale_y_continuous(labels = c(-3,0,3,6),
                     breaks = c(-0.003,0,0.003,0.006)) +
  scale_x_continuous(expand = c(0,0),
                     labels = c(-1,-0.5,0,0.5,1),
                     breaks = c(-1,-0.5,0,0.5,1),
                     limits = c(-1,1)) +
  coord_cartesian(ylim = c(-0.0035,0.01),
                  xlim = c(-1,1),
                  clip = "off")
# print(binplot)

# Fig 1c. Regression distributions

GTExRandPlot <- ggplot()
for(i in which(!CKTest$PassTest)){
  # print(i)
  deltastr <- paste0("delta",i)
  deltastr <- ensym(deltastr)
  
  GTExRandPlot <- GTExRandPlot +
    geom_line(data = randTP_GTEx,mapping = aes(x = !!deltastr,y = FST,color = "fail"),
              stat = "smooth",method = "lm",formula = y ~ poly(x,4),
              linewidth = 1,alpha = 0.1)
  
}
for(i in which(CKTest$PassTest)){
  # print(i)
  deltastr <- paste0("delta",i)
  deltastr <- ensym(deltastr)
  
  GTExRandPlot <- GTExRandPlot +
    geom_line(data = randTP_GTEx,mapping = aes(x = !!deltastr,y = FST,color = "pass"),
              stat = "smooth",method = "lm",formula = y ~ poly(x,4),
              linewidth = 1,alpha = 0.1)
  
}
GTExRandPlot <- GTExRandPlot +
  stat_smooth(data = randTP_GTEx,mapping = aes(x = deltaOG,y = FST,color = "real"),
              method = "lm",formula = y ~ poly(x,4),
              se = FALSE,linewidth = 1)

GTExRandPlot <- GTExRandPlot +
  xlab(expression("Sex difference in gene expression, \u0394")) +
  ylab(expression(atop("Between-Sex",F[ST]~"(x"*10^{-4}*")"))) +
  ggtitle("Twin Peaks pattern is consistent with absence of sex-differential selection") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.title = element_text(margin=margin(t=25,b=-25),hjust = 0.5),
        text = element_text(size = 9.5),
        plot.margin = unit(c(0, 0.15, 0, 0.15),"inches")) +
  geom_text_repel(aes(x = 0.25,y = -0.00005),color = "blue",label = "Real Data",
                  nudge_x = -0.25,nudge_y = 0.000075,size = 3.5) +
  geom_text_repel(aes(x = 0.75,y = -0.0001),color = "red",label = "Twin Peaks (21%)",
                  nudge_x = -0.5,nudge_y = -0.00005,size = 3.5) +
  geom_text_repel(aes(x = 0.75,y = -0.00015),color = "orange",segment.color = "orange",
                  label = "Non-Twin Peaks (79%)",
                  nudge_x = -0.45,nudge_y = -0.00009,size = 3.5) +
  geom_text(aes(x = -0.2,y = -0.0002),color = "black",label = "Null\nPermutations",
            size = 3.5) +
  # geom_text_repel(aes(x = 0,-0.00025),color = "black",label = "Sex-Permuted\nIterations",
  #                 nudge_x = -0.35,nudge_y = 0.00005) +
  annotate("text", x = -0.01, y = -0.00019, label = '{', angle = 0, size = 10) +
  scale_color_manual(values = c(fail = "orange",pass = "red",real = "blue"),
                     labels = c(fail = "Non-Twin Peaks Iterations",
                                pass = "Twin Peaks Iterations",
                                real = "Real Data")) +
  coord_cartesian(ylim = c(-0.0003,0.0002)) +
  scale_y_continuous(labels = c(-3,-2,-1,0,1),
                     breaks = c(-0.0003,-0.0002,-0.0001,0,0.0001)) +
  scale_x_continuous(expand = c(0,0),
                     labels = c(-1,-0.5,0,0.5,1),
                     breaks = c(-1,-0.5,0,0.5,1))
  

# print(GTExRandPlot)

fig2 <- ggarrange(GTExRandPlot,binplot,exprvarplot,
                  labels = c("A","B","C"),
                  ncol = 1,nrow = 3)

# png("./RevFigs/Fig2.png",width = 750,height = 750)
# print(fig2)
# dev.off()

# ggsave(filename = "./RevFigs/Fig2.svg",plot = fig2,
#        width = 3500,height = 3500,units = "px")

# ggsave(filename = "./RevFigs/Fig2.pdf",plot = fig2,
#        width = 3500,height = 3500,units = "px",device = cairo_pdf)
# 
# ggsave(filename = "./RevFigs/Fig2.png",plot = fig2,
#        width = 3500,height = 3500,units = "px")

ggsave(filename = "./RevFigs/Fig2_FINAL.pdf",plot = fig2,
       width = 6.5,height = 6.5,units = "in",dpi = 300,device = cairo_pdf)


# Fig 2. Weighted Regression in Tissues

pcut <- 0.05

# "liver","lung",
# "minor_salivary_gland",
# "pancreas",

tissuelist <- c("adipose_visceral_omentum","adrenal_gland",
                "artery_aorta","artery_tibial","brain_amygdala",
                "brain_caudate_basal_ganglia",
                "brain_cerebellar_hemisphere","brain_cerebellum","brain_cortex",
                "brain_frontal_cortex_ba9","brain_hippocampus",
                "brain_nucleus_accumbens_basal_ganglia",
                "brain_spinal_cord_cervical_c-1","brain_substantia_nigra","breast_mammary_tissue",
                "cells_cultured_fibroblasts","cells_ebv-transformed_lymphocytes",
                "colon_sigmoid","colon_transverse","esophagus_gastroesophageal_junction",
                "esophagus_mucosa","esophagus_muscularis","gonads","heart_atrial_appendage",
                "heart_left_ventricle","kidney_cortex","muscle_skeletal","nerve_tibial",
                "skin_not_sun_exposed_suprapubic","skin_sun_exposed_lower_leg",
                "small_intestine_terminal_ileum","stomach","thyroid","whole_blood")

tissuelistnames <- c("adipose: visceral omentum","adrenal gland",
                     "artery: aorta","artery: tibial","brain: amygdala",
                     "brain: caudate basal ganglia",
                     "brain: cerebellar hemisphere","brain: cerebellum","brain: cortex",
                     "brain: frontal cortex (BA9)","brain: hippocampus",
                     "brain: nucleus accumbens basal ganglia",
                     "brain: spinal cord cervical C-1","brain: substantia nigra","breast: mammary tissue",
                     "cells: cultured fibroblasts","cells: EBV-transformed lymphocytes",
                     "colon: sigmoid","colon: transverse","esophagus: gastroesophageal junction",
                     "esophagus: mucosa","esophagus: muscularis","gonads","heart: atrial appendage",
                     "heart: left ventricle","kidney: cortex","muscle: skeletal","nerve: tibial",
                     "skin: not sun-exposed suprapubic","skin: sun-exposed lower leg",
                     "small intestine: terminal ileum","stomach","thyroid","whole blood")

AValsDFReal <- data.frame(AVal = c(),
                          Tissue = c(),
                          pval = c(),
                          Min = c(),
                          Max = c(),
                          Name = c())

# AValsPlot <- ggplot()

for(tissue in tissuelist){
  print(tissue)
  AValsDF_i <- vroom(paste0("./AVals_New_MH/AVals_NFE.PVal",pcut,".",tissue,".MHNew.txt"),
                     show_col_types = FALSE)
  AValsDF_i$Tissue <- rep(tissue,nrow(AValsDF_i))
  
  pval <- length(which(AValsDF_i$AVal[2:1000] > AValsDF_i$AVal[1]))/1000
  
  AValsDF_i$pval <- rep(pval,1000)
  
  # AValsDFRand <- rbind(AValsDFRand,AValsDF_i[2:nrow(AValsDF_i),c("AVal","Tissue")])
  NewRow <- c(AValsDF_i[1,c("AVal","Tissue","pval")],
              quantile(AValsDF_i$AVal,0.05,names = FALSE),
              quantile(AValsDF_i$AVal,0.95,names = FALSE),
              tissuelistnames[which(tissuelist == tissue)])
  names(NewRow) <- c("AVal","Tissue","pval","Min","Max","Name")
  AValsDFReal <- rbind(AValsDFReal,NewRow)
  print(pval)
}

# nrow(AValsDFRand)
nrow(AValsDFReal)

AValsPlot <- ggplot() + geom_point(data = AValsDFReal,mapping = aes(x = AVal,y = reorder(Name,AVal))) +
  geom_vline(xintercept = 0,linetype = "dashed") +
  # geom_violin(data = AValsDFRand,mapping = aes(x = AVal,y = Tissue,color = "null"),fill = "grey",show.legend = FALSE) +
  geom_linerange(data = AValsDFReal,aes(y=Name, xmin=Min, xmax=Max,color = "null"),linewidth=1.2,show.legend = FALSE,alpha = 0.5) +
  geom_point(data = AValsDFReal,mapping = aes(x = AVal,y = Name,color = "point"),show.legend = TRUE,size = 2) +
  xlab(expression("Strength of sex-differential selection (x"*10^{-3}*"),"~italic(A))) + ylab("") +
  # ggtitle("No evidence for sex-differential selection at GTEx eQTLs") +
  scale_color_manual(values = c(null = "brown",point = "blue"),
                     labels = c(null = "Null Distribution\n90% Acceptance Region",point = "Point Estimate")) +
  # guides(color = guide_legend(title = "",
  #                             override.aes = list(shape = c(15,16),
  #                                                 color = c("brown","blue"),
  #                                                 size = c(4,2)),
  #                             byrow = TRUE)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key = element_rect(fill = "white"),legend.position = "none",
        legend.spacing.y = unit(1.5,"lines"),
        plot.title.position = "plot",
        text = element_text(size = 10)) +
  scale_x_continuous(labels = c(-2.5,0,2.5),
                     breaks = c(-0.0025,0,0.0025)) +
  geom_text_repel(aes(x = -0.0001,y = 15),color = "blue",label = "Point Estimate",
                  nudge_x = 0.0029,nudge_y = 2,size = 3.5) +
  geom_text_repel(aes(x = 0.001,y = 26),color = "brown",
                  label = "Null Distribution\n90% Acceptance Region",
                  nudge_x = 0.003,nudge_y = -3.5,size = 3.5)


# png("./RevFigs/Fig3.png",width = 700,height = 700)
# print(AValsPlot)
# dev.off()

# ggsave(filename = "./RevFigs/OLSAvals.pdf",plot = AValsPlot,
#        width = 3500,height = 3500,units = "px",device = cairo_pdf)

ggsave(filename = "./RevFigs/FigS1_FINAL.pdf",plot = AValsPlot,
       width = 6.5,height = 6.5,units = "in",dpi = 300,device = cairo_pdf)


# ggsave(filename = "./RevFigs/OLSAvals.png",plot = AValsPlot,
#        width = 3500,height = 3500,units = "px")


# Fig 3. Revision Fig

# tissuelist <- c("artery_aorta","brain_cerebellar_hemisphere","brain_cortex","brain_hippocampus",
#                 "brain_spinal_cord_cervical_c-1","cells_ebv-transformed_lymphocytes",
#                 "esophagus_mucosa","muscle_skeletal","nerve_tibial","skin_sun_exposed_lower_leg",
#                 "stomach","whole_blood")
# 
# tissuelistnames <- c("artery: aorta","brain: cerebellar hemisphere","brain: cortex",
#                      "brain: hippocampus",
#                      "brain: spinal cord cervical C-1","cells: EBV-transformed lymphocytes",
#                      "esophagus: mucosa","muscle: skeletal","nerve: tibial","skin: sun-exposed lower leg",
#                      "stomach","whole blood")


# tissuelist <- c("heart_atrial_appendage","heart_left_ventricle","esophagus_muscularis")
# tissuelist <- c("esophagus_muscularis","heart_atrial_appendage","heart_left_ventricle","esophagus_gastroesophageal_junction","colon_transverse")
# # # tissuelist <- "adipose_subcutaneous"
# tissuelistnames <- tissuelist

AValsDFRevision <- data.frame(AVal = c(),
                              Tissue = c(),
                              pval = c(),
                              Min = c(),
                              Max = c(),
                              Name = c())

# AValsPlot <- ggplot()

for(tissue in tissuelist){
  print(tissue)
  AValsDF_i <- vroom(paste0("./AVals_New_MH_Rev6/AVals_NFE.PVal",pcut,".",tissue,".MHNew.Revision.txt"),
                     show_col_types = FALSE)
  AValsDF_i$Tissue <- rep(tissue,nrow(AValsDF_i))
  
  pval <- length(which(AValsDF_i$AVal[2:1001] > AValsDF_i$AVal[1]))/1000
  
  AValsDF_i$pval <- rep(pval,1001)
  
  # AValsDFRand <- rbind(AValsDFRand,AValsDF_i[2:nrow(AValsDF_i),c("AVal","Tissue")])
  NewRow <- c(AValsDF_i[1,c("AVal","Tissue","pval")],
              quantile(AValsDF_i$AVal,0.05,names = FALSE),
              quantile(AValsDF_i$AVal,0.95,names = FALSE),
              tissuelistnames[which(tissuelist == tissue)])
  names(NewRow) <- c("AVal","Tissue","pval","Min","Max","Name")
  AValsDFRevision <- rbind(AValsDFRevision,NewRow)
  print(pval)
}

# nrow(AValsDFRand)
nrow(AValsDFRevision)

AValsPlotRev <- ggplot() + geom_point(data = AValsDFRevision,mapping = aes(x = AVal,y = reorder(Name,AVal))) +
  geom_vline(xintercept = 0,linetype = "dashed") +
  # geom_violin(data = AValsDFRand,mapping = aes(x = AVal,y = Tissue,color = "null"),fill = "grey",show.legend = FALSE) +
  geom_linerange(data = AValsDFRevision,aes(y=Name, xmin=Min, xmax=Max,color = "null"),linewidth=1.2,show.legend = FALSE,alpha = 0.5) +
  geom_point(data = AValsDFRevision,mapping = aes(x = AVal,y = Name,color = "point"),show.legend = TRUE,size = 2) +
  xlab(expression("Strength of sex-differential selection,"~italic(A))) + ylab("") +
  # ggtitle("No evidence for sex-differential selection at GTEx eQTLs") +
  scale_color_manual(values = c(null = "brown",point = "blue"),
                     labels = c(null = "Null Distribution\n90% Acceptance Region",point = "Point Estimate")) +
  # guides(color = guide_legend(title = "",
  #                             override.aes = list(shape = c(15,16),
  #                                                 color = c("brown","blue"),
  #                                                 size = c(4,2)),
  #                             byrow = TRUE)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.key = element_rect(fill = "white"),legend.position = "none",
        legend.spacing.y = unit(1.5,"lines"),
        plot.title.position = "plot",
        text = element_text(size = 10)) +
  scale_x_continuous(labels = c(-0.5,-0.25,0,0.25,0.5),
                     breaks = c(-0.5,-0.25,0,0.25,0.5)) +
  geom_text_repel(aes(x = 0,y = 16),color = "blue",label = "Point Estimate",
                  nudge_x = 0.175,nudge_y = 2,size = 3.5) +
  geom_text_repel(aes(x = 0.075,y = 27),color = "brown",
                  label = "Null Distribution\n90% Acceptance Region",
                  nudge_x = 0.15,nudge_y = -3.5,size = 3.5)

# png("TestAdiposePlot4.png",width = 700,height = 700)
# # png("presentation/NewMethodAll_wtitle.png",width = 1600,height = 1000)
# print(AValsPlotRev)
# dev.off()

# png("./RevFigs/SuppFig.png",width = 700,height = 700)
# # png("presentation/NewMethodAll_wtitle.png",width = 1600,height = 1000)
# print(AValsPlotRev)
# dev.off()

# ggsave(filename = "./RevFigs/SMAAvals.pdf",plot = AValsPlotRev,
#        width = 3500,height = 3500,units = "px",device = cairo_pdf)

ggsave(filename = "./RevFigs/Fig3_FINAL.pdf",plot = AValsPlotRev,
       width = 6.5,height = 6.5,units = "in",dpi = 300,device = cairo_pdf)


# ggsave(filename = "./RevFigs/SMAAvals.png",plot = AValsPlotRev,
#        width = 3500,height = 3500,units = "px")


CombDF <- data.frame(AVal = c(AValsDFReal$AVal[which(AValsDFReal$Name %in% AValsDFRevision$Name)],
                              AValsDFRevision$AVal[which(AValsDFRevision$Name %in% AValsDFReal$Name)]),
                    Min = c(AValsDFReal$Min[which(AValsDFReal$Name %in% AValsDFRevision$Name)],AValsDFRevision$Min[which(AValsDFRevision$Name %in% AValsDFReal$Name)]),
                    Max = c(AValsDFReal$Max[which(AValsDFReal$Name %in% AValsDFRevision$Name)],AValsDFRevision$Max[which(AValsDFRevision$Name %in% AValsDFReal$Name)]),
                    Name = c(AValsDFReal$Name[which(AValsDFReal$Name %in% AValsDFRevision$Name)],AValsDFRevision$Name[which(AValsDFRevision$Name %in% AValsDFReal$Name)]),
                    Method = c(rep("OLS",nrow(AValsDFReal[which(AValsDFReal$Name %in% AValsDFRevision$Name),])),rep("SMA",nrow(AValsDFRevision[which(AValsDFRevision$Name %in% AValsDFReal$Name),]))))

AValsPlotComb <- ggplot() + 
  geom_point(data = CombDF,
             mapping = aes(x = AVal,y = reorder(Name,AVal),shape = Method),
             position = position_dodgev(height = 1)) +
  geom_vline(xintercept = 0,linetype = "dashed") +
  # geom_violin(data = AValsDFRand,mapping = aes(x = AVal,y = Tissue,color = "null"),fill = "grey",show.legend = FALSE) +
  geom_linerange(data = CombDF,aes(y=Name, xmin=Min, xmax=Max,color = "null",linetype = Method),linewidth=1,show.legend = FALSE,alpha = 0.5,
                 position = position_dodgev(height = 1)) +
  geom_point(data = CombDF,mapping = aes(x = AVal,y = Name,color = "point",shape = Method),show.legend = TRUE,
             position = position_dodgev(height = 1)) +
  xlab(expression("Strength of sex-differential selection,"~italic(A))) + ylab("") +
  # ggtitle("No evidence for sex-differential selection at GTEx eQTLs") +
  scale_color_manual(values = c(null = "brown",point = "blue"),
                     labels = c(null = "Null Distribution\n90% Acceptance Region",point = "Point Estimate")) +
  # guides(color = guide_legend(title = "",
  #                             override.aes = list(shape = c(15,16),
  #                                                 color = c("brown","blue"),
  #                                                 size = c(4,2)),
  #                             byrow = TRUE)) +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.key = element_rect(fill = "white"),legend.position = "none",
                     legend.spacing.y = unit(1.5,"lines"),
                     plot.title.position = "plot") +
  scale_x_continuous(labels = c(-0.5,-0.25,0,0.25,0.5),
                     breaks = c(-0.5,-0.25,0,0.25,0.5)) +
  annotate("rect",ymin = seq(0.5,(nrow(CombDF) - 1)/2,by = 2),ymax = seq(1.5,(nrow(CombDF))/2,by = 2),
           xmin = -Inf,xmax = Inf,alpha = 0.1)

png("./RevFigs/CombinedRevPlot.png",width = 700,height = 700)
# png("presentation/NewMethodAll_wtitle.png",width = 1600,height = 1000)
print(AValsPlotComb)
dev.off()

CompDF <- data.frame(OldA = CombDF$AVal[which(CombDF$Method == "OLS")],
                     NewA = CombDF$AVal[which(CombDF$Method == "SMA")],
                     OldMin = CombDF$Min[which(CombDF$Method == "OLS")],
                     OldMax = CombDF$Max[which(CombDF$Method == "OLS")],
                     NewMin = CombDF$Min[which(CombDF$Method == "SMA")],
                     NewMax = CombDF$Max[which(CombDF$Method == "SMA")])

ComparisonFig <- ggplot() +
  geom_linerange(data = CompDF,aes(y=NewA, xmin=OldMin, xmax=OldMax),
                 color = "grey") +
  geom_linerange(data = CompDF,aes(x=OldA, ymin=NewMin, ymax=NewMax),
                 color = "grey") +
  geom_point(data = CompDF,
             mapping = aes(x = OldA,y = NewA)) +
  geom_abline(intercept = 0,slope = 1,linetype = "dashed") +
  theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.key = element_rect(fill = "white"),legend.position = "none",
                     legend.spacing.y = unit(1.5,"lines"),
                     plot.title.position = "plot") +
  xlab("OLS (original method) AVals") + ylab("SMA (new method) AVals")

png("./RevFigs/ComparisonPlot.png",width = 500,height = 300)
# png("presentation/NewMethodAll_wtitle.png",width = 1600,height = 1000)
print(ComparisonFig)
dev.off()

