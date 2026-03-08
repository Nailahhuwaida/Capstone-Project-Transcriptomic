# Capstone-Project-Transcriptomic (RSCRIPT)
#Analisis GEN ROS-related Sebagai Potensi Biomarker Diagnostik Pada Tuberculosis
#Dataset: GSE313408 (Tuberculosis)
#Platform: Microarray (Affymetrix GPL31250)
#Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG) 
#Differential expression analysis with limma

if (!require("BiocManager", quietly = TRUE))  {
  install.packages("BiocManager") 
}

BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE) 

install.packages(c("pheatmap", "ggplot2", "dplyr"))

#umap: grafik (plot UMAP) 
if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}

#heatmap
install.packages("pheatmap")

#Memanggil library 
#library() digunakan agar fungsi di dalam package bisa digunakan 
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(umap)

# load series and platform data from GEO

gset <- getGEO("GSE313408", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL31250", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- "11111111110000000000"
sml <- strsplit(gsms, split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("Healthy","Disease"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ] # skip missing values

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","SEQUENCE","SPOT_ID"))
write.table(tT, file=stdout(), row.names=F, sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# summarize test results as "up", "down" or "not expressed"
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.05, lfc=0)

# Venn diagram of results
vennDiagram(dT, circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good], fit2$df.total[t.good], main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
# Please note that the code provided to generate graphs serves as a guidance to
# the users. It does not replicate the exact GEO2R web display due to multitude
# of graphical options.
# 
# The following will produce basic volcano plot using limma function:
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2, column=ct, status=dT[,ct], legend=F, pch=20, cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
ord <- order(gs)  # order samples by group
palette(c("#1B9E77", "#7570B3", "#E7298A", "#E6AB02", "#D95F02",
          "#66A61E", "#A6761D", "#B32424", "#B324B3", "#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE313408", "/", annotation(gset), sep ="")
boxplot(ex[,ord], boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=gs[ord])
legend("topright", groups, fill=palette(), bty="n")

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE313408", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, group=gs, main=title, legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 9, random_state = 123)
par(mar=c(3,3,2,6), xpd=TRUE)
plot(ump$layout, main="UMAP plot, nbrs=9", xlab="", ylab="", col=gs, pch=20, cex=1.5)
legend("topright", inset=c(-0.15,0), legend=levels(gs), pch=20,
       col=1:nlevels(gs), title="Group", pt.cex=1.5)


# mean-variance trend, helps to see if precision weights are needed
plotSA(fit2, main="Mean variance trend, GSE313408")

#VISUALISASI HEATMAP 
topTableResults <- tT[order(tT$adj.P.Val), ]

top50 <- head(topTableResults, 50)

mat_heatmap <- ex[top50$ID, ]

rownames(mat_heatmap) <- top50$ID

mat_heatmap <- mat_heatmap[rowSums(is.na(mat_heatmap)) == 0, ]

gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]

annotation_col <- data.frame(Group = gset$group)
rownames(annotation_col) <- colnames(mat_heatmap)

pheatmap(
  mat_heatmap,
  scale = "row",                 # Z-score per gen
  annotation_col = annotation_col,
  show_colnames = FALSE,         # nama sampel dimatikan
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)

#Analisis Enrichment

# 1. Install & Load Library Tambahan
if (!require("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!require("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!require("enrichplot", quietly = TRUE)) BiocManager::install("enrichplot")

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

annotation <- fData(gset)

topTableResults$GENE_ID <- annotation$GENE_ID[
  match(topTableResults$ID, annotation$ID)
]

topTableResults$GENE_SYMBOL <- annotation$GENE_SYMBOL[
  match(topTableResults$ID, annotation$ID)
]

deg_signifikan <- topTableResults[
  abs(topTableResults$logFC) > 1 & topTableResults$adj.P.Val < 0.05,
]

genes_to_test <- deg_signifikan$GENE_ID
genes_to_test <- genes_to_test[!is.na(genes_to_test)]
genes_to_test <- unique(genes_to_test)

ego_all <- enrichGO(
  gene          = genes_to_test,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

dotplot(ego_all, showCategory = 10, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free") +
  ggtitle("Gene Ontology Enrichment (BP, CC, MF)")


# KEGG Pathway
kegg_res <- enrichKEGG(
  gene = genes_to_test,
  organism = "hsa",
  pvalueCutoff = 0.1
)

barplot(kegg_res, showCategory = 10) +
  ggtitle("KEGG Pathway Enrichment")


#MENYIMPAN HASIL 

# write.csv(): menyimpan hasil analisis ke file CSV
write.csv(topTableResults, "Hasil_GSE313408_DEG.csv")

message("Analisis selesai. File hasil telah disimpan.")


