# Cheng Rep

itnum <- commandArgs(trailingOnly = TRUE)[1]

norm.counts.filt <- read.table("./norm_counts_filt.Ensembl77.GTEx_v3.txt",header = TRUE,
                               stringsAsFactors = FALSE)

randTP <- read.table("./AncTP77v3RANDAll/AllAncTP.RAND_GTEx.avgHUD.All.txt",
                     header = TRUE,stringsAsFactors = FALSE)

norm.counts.filt <- norm.counts.filt[which(norm.counts.filt$GeneName %in% randTP$name),]

norm.counts.filt[,3:22] <- norm.counts.filt[,sample(3:22)]
norm.counts.filt$females <- rowMeans(norm.counts.filt[,3:8])
norm.counts.filt$males <- rowMeans(norm.counts.filt[,9:22])
norm.counts.filt$delta <- (norm.counts.filt$males - norm.counts.filt$females)/(norm.counts.filt$males + norm.counts.filt$females)

new1 <- paste0("delta",itnum)

randTP[[new1]] <- norm.counts.filt$delta

write.table(randTP,"./AncTP77v3RANDAll/AllAncTP.RAND_GTEx.avgHUD.All.txt",
            col.names = TRUE,row.names = FALSE,quote = FALSE)
