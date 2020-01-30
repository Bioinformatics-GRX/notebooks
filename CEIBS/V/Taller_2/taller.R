library(ggbiplot)
library(factoextra)
library(cluster)
library(fpc)
library(gplots)
library(pheatmap)
library(edgeR)
library(topGO)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# 0. Load Input File 
geneExpressionMatrix_file <- "SKCM_all_normalized.normalized_results.tsv"
geneExpressionMatrix <- read.delim(geneExpressionMatrix_file,stringsAsFactors = F)
geneExpressionMatrix <- geneExpressionMatrix[grep("?",geneExpressionMatrix$gene_id,fixed = T,invert = T),]
rownames(geneExpressionMatrix) <- geneExpressionMatrix$gene_id
geneExpressionMatrix$gene_id <- NULL
colnames(geneExpressionMatrix) <- paste0("P",c(1:ncol(geneExpressionMatrix)))

# 1. Filter and Subset 
geneVariance <- apply(geneExpressionMatrix,1,var)
patientVariance <- apply(geneExpressionMatrix,2,var)

geneMOSTvariants <- which(geneVariance > mean(geneVariance))
patientsMOSTvariants <- which(patientVariance > mean(patientVariance))
filteredMatrix <- geneExpressionMatrix[geneMOSTvariants,patientsMOSTvariants]
rownames(filteredMatrix) <- gsub("\\|.*","",rownames(filteredMatrix))

# Some HeatMaps
#filteredMatrix <- na.omit(scale(filteredMatrix))
scpf <- pheatmap(filteredMatrix,scale = "row")
save_pheatmap_pdf(scpf, "filteredMatrix.pdf",55,20)

# 2. Clutering of Patients
  # 2.1. How many clusters should we expect?
wss <- (nrow(filteredMatrix)-1)*sum(apply(filteredMatrix,2,var))
for (i in 2:15) {
  wss[i] <- sum(kmeans(filteredMatrix,centers=i)$withinss)
}
plot(1:15, wss, type="b", xlab="Number of Clusters",ylab="Within groups sum of squares")

  # 2.2. K-means 
km_fit_patients <- kmeans(filteredMatrix, 4) # Centroides
fviz_cluster(km_fit_patients, data = filteredMatrix)

patientsASrows <- t(filteredMatrix)
km_fit_patients <- kmeans(patientsASrows, 4) # Centroides
fviz_cluster(km_fit_patients, data = patientsASrows)

outlayer <- table(km_fit_patients$cluster) == 1
outlayer <- names(outlayer)[outlayer] 
outlayer <- which(km_fit_patients$cluster == outlayer)
patientsASrows <- patientsASrows[-outlayer,]

km_fit_patients <- kmeans(patientsASrows, 3) # Centroides
fviz_cluster(km_fit_patients, data = patientsASrows)
table(km_fit_patients$cluster)
km_fit_patients$cluster[1:5]

## 3. Differential Expression among Groups of Patients
DEobj <- DGEList(t(patientsASrows))
groups <- as.factor(km_fit_patients$cluster)
plotMDS(DEobj, col = as.numeric(groups))
mm <- model.matrix(~0 + groups)
y <- voom(t(patientsASrows), mm)

fit <- lmFit(y, mm)
head(coef(fit))

contr12 <- makeContrasts(groups1 - groups2, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr12)
tmp <- eBayes(tmp)
top.table12 <- topTable(tmp, sort.by = "P", n = Inf)

contr13 <- makeContrasts(groups1 - groups3, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr13)
tmp <- eBayes(tmp)
top.table13 <- topTable(tmp, sort.by = "P", n = Inf)

contr23 <- makeContrasts(groups2 - groups3, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr23)
tmp <- eBayes(tmp)
top.table23 <- topTable(tmp, sort.by = "P", n = Inf)

commontopDiffGenes <- unique(c(row.names(top.table12),row.names(top.table13),row.names(top.table23)))
pheatmap(filteredMatrix[commontopDiffGenes,],scale = "row", annotation_col = data.frame(groups))

# 4. Functional Annotation
geneUniverse <- row.names(filteredMatrix)
# 4.1. Group 1 vs 2
genes12names <- gsub("\\|.*","",row.names(top.table12))
matching <- match(geneUniverse,genes12names) 
geneUniverse12 <- top.table12$adj.P.Val[matching]
names(geneUniverse12) <- geneUniverse

GOdata <- new("topGOdata", ontology = "BP",
           allGenes = geneUniverse12, geneSel = topDiffGenes,
           nodeSize = 1, mapping  = "org.Hs.eg.db",
           annot = annFUN.org, ID="symbol")

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

allRes12 <- GenTable(GOdata, classicFisher = resultFisher, classicKS = resultKS,
                    elimKS = resultKS.elim, orderBy = "elimKS", ranksOf = "classicFisher", 
                    topNodes = 25)

showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')
# 4.2. Group 1 vs 2
genes13names <- gsub("\\|.*","",row.names(top.table13))
matching <- match(geneUniverse,genes13names) 
geneUniverse13 <- top.table13$adj.P.Val[matching]
names(geneUniverse13) <- geneUniverse

GOdata <- new("topGOdata", ontology = "BP",
              allGenes = geneUniverse13, geneSel = topDiffGenes,
              nodeSize = 1, mapping  = "org.Hs.eg.db",
              annot = annFUN.org, ID="symbol")

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

allRes13 <- GenTable(GOdata, classicFisher = resultFisher, classicKS = resultKS,
                     elimKS = resultKS.elim, orderBy = "elimKS", ranksOf = "classicFisher", 
                     topNodes = 25)

showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')

# 4.2. Group 2 vs 3
genes23names <- gsub("\\|.*","",row.names(top.table13))
matching <- match(geneUniverse,genes23names) 
geneUniverse23 <- top.table23$adj.P.Val[matching]
names(geneUniverse23) <- geneUniverse

GOdata <- new("topGOdata", ontology = "BP",
              allGenes = geneUniverse23, geneSel = topDiffGenes,
              nodeSize = 1, mapping  = "org.Hs.eg.db",
              annot = annFUN.org, ID="symbol")

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

allRes23 <- GenTable(GOdata, classicFisher = resultFisher, classicKS = resultKS,
                     elimKS = resultKS.elim, orderBy = "elimKS", ranksOf = "classicFisher", 
                     topNodes = 25)

showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')
