################################################################################
# Cervical Cancer vs Normal Gene Expression Analysis
# CGC Congress 2023 - GSE63514 Dataset
#
# Description: This script performs the following tasks:
#         1) Load and normalize Affymetrix microarray data
#         2) Sample selection for Normal and Cancer groups
#         3) Quality control and PCA visualization
#         4) Differential expression analysis using limma
#         5) DEG visualization: Volcano plot and heatmap
# 
# Note: CEL files must be present in the specified directory before running the script.
################################################################################

#--------------------- Load Libraries ---------------------
library(affy)
library(limma)
library(annotate) 
library(hgu133plus2.db) 
library(pheatmap)
library(gplots)
library(org.Hs.eg.db)
library(ggplot2)

#--------------------- Load Data ---------------------
# Set working directory for CEL files
setwd("D:/Project/Cervical/datac") # Change to the directory containing CEL files
# Read Affymetrix CEL files into a gset object
gset <- ReadAffy() 
pData(gset) 

#--------------------- Normalization ---------------------
### Normalization using RMA (Robust Multi-array Average)
data <- rma(gset)  
ex <- exprs(data)  
# Check for log2 transformation needs
max(ex)  
min(ex)  
# ex <- log2(ex + 1)  Commented out the log2 transformation as it was not needed for this project

#--------------------- Sample Selection ---------------------
setwd("D:/Project/Cervical/result") 
# Define sample groups: X = Excluded, 1 = Normal, 0 = Cancer
gsms <- paste0("111111111111111111111111XXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "0000000000000000000000000000")
sml <- strsplit(gsms, split = "")[[1]]  
sel <- which(sml != "X")  
sml <- sml[sel]  
data <- data[, sel] 
ex <- exprs(data)  
dim(ex) 

# Define group factor: 24 Normal and 28 Cancer
gr <- factor(c(rep("Normal", 24), rep("Cancer", 28)))  
#Save Data:
setwd("D:/Project/Cervical/RData")
save(data, ex, gr, file = "EXDATA.RData")

#--------------------- QC Plots ---------------------
setwd("D:/Project/Cervical/plots")  

# Boxplot for expression data
png("boxplot_N_vs_C.png")
boxplot(ex, main = "Expression Data Boxplot")
dev.off()

# PCA (Principal Component Analysis)
# Method 1: PCA - Row means centering
exm <- ex - rowMeans(ex)  
pc <- prcomp(exm)  
pcr <- data.frame(pc$r[, 1:3], gr)  
png("PCA_rowMeans_N_vs_C.png", width = 2000, height = 1500, res = 300, units = "px")  
ggplot(pcr, aes(PC1, PC2, color = gr)) + ggtitle("PCA rowMeans") +
  geom_label(size = 2.2, fontface = "bold", label = gr) +  
  geom_vline(xintercept = 0, col = "black", linetype = "longdash") +  
  geom_hline(yintercept = 0, col = "black", linetype = "longdash") + 
  stat_ellipse(aes(fill = gr), type = "norm", geom = "polygon", alpha = 0.09, level = 0.79) +  
  theme_bw()  
dev.off() 

# Method 2: PCA using scaled data
ex.scale <- t(scale(t((ex)), scale = TRUE, center = FALSE))  
pc <- prcomp(ex.scale) 
pcr <- data.frame(pc$r[, 1:3], Group = gr)  
png("PCA_scale_N_vs_C.png", width = 2000, height = 1500, res = 300, units = "px")  
ggplot(pcr, aes(PC1, PC2, color = gr)) + ggtitle("PCA scaled") +
  geom_label(size = 2.2, fontface = "bold", label = gr) +  
  geom_vline(xintercept = 0, col = "black", linetype = "longdash") +  
  geom_hline(yintercept = 0, col = "black", linetype = "longdash") + 
  stat_ellipse(aes(fill = gr), type = "norm", geom = "polygon", alpha = 0.09, level = 0.79) +  
  theme_bw() 
dev.off()  

#--------------------- Differential Expression Analysis ---------------------
setwd("D:/Project/Cervical/result")

# Design matrix for differential expression
gs <- factor(sml)  # X = Excluded, 0 = Cancer, 1 = Normal
groups <- make.names(c("Cancer","Normal"))  
levels(gs) <- groups
data$group <- gs  
design <- model.matrix(~group + 0, data)
colnames(design) <- levels(gs)  

# Linear modeling with limma
fit <- lmFit(data, design)
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust.method="fdr", sort.by="logFC", number=Inf)

# Annotate DEGs with gene symbols
probeList <- rownames(tT)
Gene.symbol <- getSYMBOL(probeList, 'hgu133plus2.db')
top <- cbind(tT, Gene.symbol)
top <- na.omit(top)

# Filter significant DEGs
degs <- subset(top, abs(logFC) > 1.5 & adj.P.Val < 0.05)
degs$ID = rownames(degs) 

# Aggregate DEGs by Gene.symbol
nums <- table(degs[,"Gene.symbol"])  
degs.ag <- data.frame(matrix(0, nrow=length(nums), ncol=ncol(degs)))  
for(i in 1:length(nums)) {
  j <- which(degs[,"Gene.symbol"] == names(nums[i]))
  fac <- degs[j, "logFC"]  
  k <- which(max(fac) == fac) 
  degs.ag[i,] <- degs[j[k],] 
}
colnames(degs.ag) <- colnames(degs)
rownames(degs.ag) <- degs.ag[,"Gene.symbol"]

# Count upregulated and downregulated DEGs
degsup <- subset(degs.ag, logFC > 1.5 & adj.P.Val < 0.05)  # Upregulated
degsdown <- subset(degs.ag, logFC < -1.5 & adj.P.Val < 0.05)  # Downregulated

# Save aggregated DEGs for downstream analysis
write.table(degs.ag, file="Degs - aggregated in R.txt", quote=F, sep="\t")
# Write DEGs' gene symbols for STRING database input
write.table(degs.ag$Gene.symbol, file="Degs GeneNames.txt", quote=F, row.names=F, col.names=F)
# Save data objects for later use
setwd("D:/Project/Cervical/RData")
save(degs, degs.ag, degsdown, degsup, design, file="DEGs.RData")

#--------------------- Volcano Plot for DEGs ---------------------
setwd("D:/Project/Cervical/plots")
volc <- subset(top)
volc$Significant <- "No"
volc$Significant[volc$logFC > 1.5 & volc$adj.P.Val < 0.05] <- "Up"
volc$Significant[volc$logFC < -1.5 & volc$adj.P.Val < 0.05] <- "Down"

# Create volcano plot
png("volcano plot.png", height=1600, width=2000, res=300, units="px")
ggplot(volc, aes(logFC, -log10(adj.P.Val), color=Significant)) +
  geom_point(size=1.6, shape=19) + theme_bw() +
  geom_vline(xintercept=c(-1.5, 1.5), col="#EF00FF", linetype="longdash") +
  geom_hline(yintercept=-log10(0.05), col="#EF00FF", linetype="longdash") +
  scale_color_manual(values=c("blue", "gray", "red")) +
  ggtitle("Volcano Plot")
dev.off()

#--------------------- Heatmap for DEGs ---------------------
# Convert expression matrix to a dataframe
bayan <- data.frame(ex) 
heatmap_data <- bayan[rownames(bayan) %in% degs.ag$ID,]  
#making annotation_col:
coloo <- colnames(heatmap_data) 
annotation1 <- data.frame(gr, coloo)
rownames(annotation1) <- annotation1$coloo
colnames(annotation1) = "Group"
annotation1 <- annotation1[,-2,drop=F]
annotation2 <- data.frame(degs.ag$logFC)
rownames(annotation2) <- degs.ag$Gene.symbol
colnames(annotation2) = "Up/Down regulated"
annotation2$`Up/Down regulated`[degs.ag$logFC > 1.5 & degs.ag$adj.P.Val < 0.05] <- "Up"
annotation2$`Up/Down regulated`[degs.ag$logFC < -1.5 & degs.ag$adj.P.Val < 0.05] <- "Down"
ann_colors = list(
  "Group" = c( Normal= "#1B9E77",  Cancer= "#D95F02"),
  "Up/Down regulated" = c(Down = "#7570B3", Up = "#E7298A"))

# Use pheatmap for heatmap generation
setwd("D:/Project/Cervical/plots")
png("heatmap_DEGs.png",height = 3500, width = 3300,res=300,units = "px" )
pheatmap(heatmap_data,color= greenred(256), show_rownames=F ,border_color =NA
         ,annotation_col = annotation1,annotation_row = annotation2, 
         scale = "row", labels_col = gr, annotation_colors = ann_colors,
         annotation_names_row = F,annotation_names_col = F,fontsize_row = 4.5,
         fontsize_col=5,angle_col=45
         clustering_distance_rows = "correlation",
         clustering_distance_cols= "correlation",
         clustering_method = "complete")
dev.off()

################################################################################
# End of Script
# Summary: This analysis identified and visualized differentially expressed genes
#          between Cervical Cancer and Normal samples, with results saved for 
#          downstream analysis and visualization.
################################################################################

