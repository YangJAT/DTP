library(Seurat)
library(plyr)
library(future)
library(ggplot2)
library(cowplot)
library(data.table)
library(tidyverse)
library(future.apply)
library(RColorBrewer)

setwd("/home/ug0302/workspace/CMN")
source("function.R")

path = "/home/ug0302/workspace/CMN/rawdata"
datafilt <- read10x_all(path = path, merge = T)

name = "QC_plot.pdf"
pdf(file = name, width = 10, height = 7)
qc_plot(data = datafilt)
dev.off()

datafilt <- qc_filter(datafilt, 
                      nFeature_dn = 1000, nFeature_up = 10000,
                      mito_thres = 20, gene_filter = 3)

datafilt <- doublet_rm(datafilt, prop_doublet = 0.05)
plot <- dimplot_new(datafilt, pt.size = 1, label = T, group.by = c("DFclass"))
ggsave("umap_doublet.png", plot, dpi = 300, width = 6, height = 5)
datafilt <- datafilt[,datafilt$DFclass == "Singlet"]

load("/home/ug0302/CITEseq/public_data/gtf_v22.Rdata")
common <- intersect(rownames(datafilt), gene_info$SYMBOL)
datafilt <- datafilt[common,]

datafilt <- autocluster(datafilt)

plot <- dimplot_new(datafilt, pt.size = 1, label = T,
                    group.by = c("seurat_clusters"))
ggsave("umap_cluster_all.png", plot, dpi = 300, width = 6, height = 5)

plot <- dimplot_new(datafilt, pt.size = 1, label = T,
                    group.by = c("sample"))
ggsave("umap_sample_all.png", plot, dpi = 300, width = 6.5, height = 5)


# CytoTRACE

library(CytoTRACE)
library(tibble)

input <- GetAssayData(datafilt, slot = "counts",
                      assay = "RNA") %>% as.matrix()

score <- CytoTRACE(input, enableFast = FALSE)
score <- data.frame(id = names(score[["CytoTRACE"]]),
                    cytotrace = unname(score[["CytoTRACE"]]))

score <- column_to_rownames(score, var = "id")
datafilt <- AddMetaData(datafilt, score)

# Cell entropy

library(vegan)
library(tibble)

input <- GetAssayData(datafilt, slot = "counts",
                      assay = "RNA") %>% as.matrix()

# entropy

score <- do.call(rbind, lapply(colnames(input), function(i){
  data <- input[,i]
  equ <- diversity(data)
  data.frame(id = i, entropy = equ)
}))

score <- column_to_rownames(score, var = "id")
datafilt <- AddMetaData(datafilt, score)


# mRNAsi 

data <- GetAssayData(datafilt, slot = "counts",
                     assay = "RNA") %>% as.matrix()

path = "/home/ug0302/CITEseq/public_data/pcbc-stemsig.tsv"
coef <- read.delim(path, header = FALSE, row.names = 1 ) %>%
  as.matrix() %>% drop()

common <- intersect(rownames(data), names(coef))
data <- data[common,]
coef <- coef[rownames(data)]

score <- apply(data, 2, function(z){cor(z, coef,
                                        method="sp", use="complete.obs")})
score <- score - min(score)
score <- score / max(score)
score <- data.frame(mRNAsi = score)
datafilt <- AddMetaData(datafilt, score)


input <- data.frame(type = subfilt$seurat_clusters,
                    value = subfilt$prolif_score)

name = "boxplot.pdf"
common_boxplot2(input = input,
                output = name,
                show.point = TRUE, order = TRUE,
                decreasing = TRUE, rotate = 45,
                width = 7, height = 5)



# cell cycle 

cellcycle <- cellcycle_seurat(datafilt)
datafilt <- AddMetaData(datafilt, cellcycle)
cor.test(datafilt$mRNAsi, datafilt$prolif_score)

plot <- dimplot_new(datafilt, pt.size = 1, label = T,
                    group.by = c("cyc_tri"))
ggsave("umap_cycle_all.png", plot, dpi = 300, width = 7, height = 5)



path = "/home/ug0302/workspace/CMN/tumor_sig1.csv"
scores <- seurat_score(datafilt, source = path,
                       geneset = NULL,
                       min.sz = 10)
datafilt <- AddMetaData(datafilt, scores)


select_al <- colnames(scores)
all_plots <- featureplot_new(datafilt, pt.size = 1.5,
                             reduction = "umap", color = "blue2red",
                             features = select_al, outlier.rm = F)

name = "umap_features.png"
export_featureplot(all_plots = all_plots,
                   ncol = 3, dpi = 300, output = name)


subfilt <- datafilt[,datafilt$sample == "PS002n"]
subfilt <- autocluster(subfilt, nfeatures = 1000, ndim = 15,
                       neigh = 20, dist = 0.5, res = 0.3)

plot <- dimplot_new(subfilt, pt.size = 1, label = T,
                    group.by = c("seurat_clusters"))
ggsave("umap_cluster_002.png", plot, dpi = 300, width = 6, height = 5)

plot <- dimplot_new(subfilt, pt.size = 1, label = T,
                    group.by = c("cyc_tri"))
ggsave("umap_cycle_002.png", plot, dpi = 300, width = 7, height = 5)


select_al <- c("cytotrace", "entropy", "mRNAsi")

all_plots <- featureplot_new(subfilt, pt.size = 1.5,
                             reduction = "umap", color = "blue2red",
                             features = select_al, outlier.rm = F)

name = "stem_within002.png"
export_featureplot(all_plots = all_plots,
                   ncol = 3, dpi = 300, output = name)

diff <- seurat_diffall(subfilt, min.pct = 0.25,
                       group.by = "cyc_tri")

diff_flt <- diff[(diff$logfc > 0.5 & diff$adjp < 0.05),]
write.table(diff_flt, "diff_within002.txt", sep = "\t",
            quote = F, row.names = F, col.names = T)

input <- data.frame(gene = diff_flt$gene,
                    type = diff_flt$cluster)

path = "/home/ug0302/CITEseq/public_data/all_genesets.rds"
enrich <- enrichr_cluster(data = input,
                          source = path, geneset = "msigdb_hallmark",
                          cutoff_p = 0.05, cutoff_adjustp = 1)

write.table(enrich, "enrich_within002.txt", sep = "\t",
            quote = F, row.names = F, col.names = T)

name = "heatmap_enrich_within002.pdf"
pdf(file = name, width = 6, height = 7)
enrichr_cluster_plot(enrich, topN = 5, 
                     select = "pvalue", color = "col2")
dev.off()


path = "/home/ug0302/workspace/CMN/tumor_sig1.csv"
scores <- seurat_score(subfilt, source = path,
                       geneset = NULL,
                       min.sz = 10)
subfilt <- AddMetaData(subfilt, scores)


select_al <- colnames(scores)
all_plots <- featureplot_new(subfilt, pt.size = 1.5,
                             reduction = "umap", color = "blue2red",
                             features = select_al, outlier.rm = F)

name = "umap_features_within002.png"
export_featureplot(all_plots = all_plots,
                   ncol = 3, dpi = 300, output = name)



# sample-wise ==================================================================

unique(datafilt$sample)
datafilt$sample[datafilt$sample == "PS003n" |
                  datafilt$sample == "PS005n"] <- "PS003_5n"
datafilt$sample[datafilt$sample == "PS004n" |
                  datafilt$sample == "PS006n"] <- "PS004_6n"

diff <- seurat_diffall(datafilt, min.pct = 0.25,
                       group.by = "sample")

diff_flt <- diff[(diff$logfc > 0.5 & diff$adjp < 0.05),]
write.table(diff_flt, "diff_sample.txt", sep = "\t",
            quote = F, row.names = F, col.names = T)

# enrichment

input <- data.frame(gene = diff_flt$gene,
                    type = diff_flt$cluster)

path = "/home/ug0302/CITEseq/public_data/all_genesets.rds"
enrich <- enrichr_cluster(data = input,
                          source = path, geneset = "KEGG",
                          cutoff_p = 0.05, cutoff_adjustp = 1)

write.table(enrich, "enrich_sample.txt", sep = "\t",
            quote = F, row.names = F, col.names = T)

name = "heatmap_enrich_sample.pdf"
pdf(file = name, width = 6, height = 7)
enrichr_cluster_plot(enrich, topN = 5, 
                     select = "pvalue", color = "col2")
dev.off()





