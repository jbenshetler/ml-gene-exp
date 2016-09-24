# http://bioinformatics.knowledgeblog.org/2011/06/20/analysing-microarray-data-in-bioconductor/
biocLite("hgu133plus2.db")
install.packages("ggplot2")


setwd("C:/R/geo/data/GSE32269")
library(simpleaffy)
library(RColorBrewer)
library(affyPLM)
library(hgu133plus2.db)
library(annotate)
library(ggplot2)

celfiles <- read.affy(covdesc="phenodatashort.txt", path="C:/R/geo/data/GSE32269")
celfiles
celfiles.gcrma <- gcrma(celfiles)

cols <- brewer.pal(8, "Set1")
boxplot(celfiles, col=cols)
boxplot(celfiles.gcrma, col=cols)
hist(celfiles, col=cols)
hist(celfiles.gcrma, col=cols)
celfiles.qc <- fitPLM(celfiles)
image(celfiles.qc, which=1, add.legend=TRUE)
image(celfiles.qc, which=4, add.legend=TRUE)
RLE(celfiles.qc, main="RLE")
NUSE(celfiles.qc, main="NUSE")
NUSE(celfiles.qc, main="NUSE",las=2)
eset <- exprs(celfiles.gcrma)
distance <- dist(t(eset),method="maximum")
clusters <- hclust(distance)
plot(clusters)

# Filter out uninteresting 
celfiles.filtered <- nsFilter(celfiles.gcrma, require.entrez=FALSE, remove.dupEntrez=FALSE)

samples <- celfiles.gcrma$Target
samples <- as.factor(samples)
design <- model.matrix(~0 + samples)
colnames(design) <- c("normal", "cancer")

library(limma)
fit <- lmFit(exprs(celfiles.filtered$eset), design)

contrast.matrix <- makeContrasts(cancer_normal = cancer - normal, levels=design)
cancer_fits <- contrasts.fit(fit, contrast.matrix)
cancer_ebfit <- eBayes(cancer_fits)

topTable(cancer_ebfit, number=20, coef=1)

nrow(topTable(cancer_ebfit, coef=1, number=10000, lfc=5))

gene_list <- topTable(cancer_ebfit, coef=1, number=1000000, sort.by="logFC")
head(gene_list$logFC)
head(gene_list$P.Value)

par(mar=c(3,3,2,1), mgp=c(2,.7,0), tck=-.01)

plot(gene_list$logFC, -log10(gene_list$P.Value),
     xlim=c(-10, 10), ylim=c(0, 15), #Set limits
     xlab="log2 fold change", ylab="-log10 p-value")#Set axis labels


sum(abs(gene_list$logFC) > 2 & gene_list$P.Value < 0.001)

no_of_genes = dim(probeset.list)[1]
sum(abs(gene_list$logFC) > 2 & gene_list$P.Value < 0.05/no_of_genes)

sum(abs(gene_list$logFC) > 2 & gene_list$adj.P.Val < 0.05)

#gene_list$threshold = as.factor(abs(gene_list$logFC) > 2 & gene_list$P.Value < 0.05/no_of_genes)

gene_list$threshold = as.factor(abs(gene_list$logFC) > 2 & gene_list$P.Value < 0.05)

g = ggplot(data=gene_list, aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-10, 10)) + ylim(c(0, 15)) +
  xlab("log2 fold change") + ylab("-log10 p-value")

g
