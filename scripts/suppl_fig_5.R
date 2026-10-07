## Supplementary 5
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(RColorBrewer)


TSK <- readLines('data/TSK_markers.txt')
patients <- c('4','6','16')
cutoff <- 0.5
Seurat_counts_list <- list()
for (Patient in patients){
  pattern <- paste0("ST_P",Patient,"_rep(\\d+)")
  Samples <- dir('data/VISIUM/',pattern = pattern)[sapply(dir('data/VISIUM/',pattern = pattern), function(x) dir.exists(paste0('data/VISIUM/',x,'/spatial')))]
  for (Sample in Samples){
    skip_post <- F
    print(Sample)
    patient_folder <- paste0('data/VISIUM/',Sample)
    Seurat_counts <- Load10X_Spatial(
      paste0(patient_folder),
      slice = Sample,
      image = Read10X_Image(paste0(patient_folder, "/spatial"))
    )
    Seurat_counts@images[[1]]@scale.factors[["lowres"]] <-  ScaleFactors(Seurat_counts@images[[1]])$hires
    ## link barcodes to prediction values (probs)
    probs <- read.table(paste0(patient_folder,'/TPM_GNN_probabilities_scaled.tsv'), header = T)
    
    poslist <- read.csv(paste0(patient_folder,'/spatial/tissue_positions_list.csv'), header = F, row.names = 1)
    colnames(poslist)[c(2,3)]<- c('x','y')
    poslist$barcodes <- rownames(poslist)
    probs <- merge(probs, poslist, by = c('x','y'))
    probs <- probs[order(match(probs$barcodes,Cells(Seurat_counts))),]
    Seurat_counts <- AddMetaData(Seurat_counts, probs$prob, col.name = 'prob')
    Seurat_counts <- NormalizeData(Seurat_counts,normalization.method = 'LogNormalize', verbose = F,scale.factor = 1e6)
    Seurat_counts <- ScaleData(Seurat_counts, verbose = F)
    Seurat_counts <- FindVariableFeatures(Seurat_counts, verbose = F)
    Seurat_counts <- RunPCA(Seurat_counts, verbose = FALSE)
    Seurat_counts <- FindNeighbors(Seurat_counts, reduction = "pca", dims = 1:30, verbose = F)
    Seurat_counts <- FindClusters(Seurat_counts, verbose = FALSE)
    Seurat_counts <- RunUMAP(Seurat_counts, reduction = "pca", dims = 1:30, verbose = F)
    Seurat_counts_list[[Sample]] <- AddModuleScore(Seurat_counts, list(TSK), name = 'TSK')
    
  }
}
cor.test(Reduce(c,sapply(Seurat_counts_list, FUN = function(x) x$prob)),Reduce(c,sapply(Seurat_counts_list, FUN = function(x) x$TSK1)), method = 'spearman')
for (samp in names(Seurat_counts_list)){
  print(samp)
  print(cor.test(Seurat_counts_list[[samp]]$prob, Seurat_counts_list[[samp]]$TSK1, method = 'spearman'))
}

translation <- c(ST_P4_rep1 = 'cSCC1',
                 ST_P4_rep2 = 'cSCC2',
                 ST_P6_rep1 = 'cSCC3',
                 ST_P6_rep2 = 'cSCC4',
                 ST_P16_rep1 = 'cSCC5',
                 ST_P16_rep2 = 'cSCC6',
                 ST_P16_rep3 = 'cSCC7',
                 ST_P16_rep4 = 'cSCC8')
# Build long dataframe
df <- bind_rows(
  lapply(names(Seurat_counts_list), function(nm) {
    data.frame(
      prob = Seurat_counts_list[[nm]]$prob,
      TSK1 = Seurat_counts_list[[nm]]$TSK1,
      sample = translation[nm]
    )
  })
)

# Add combined facet
df_combined <- df %>%
  mutate(sample = as.character(sample))

df_all <- bind_rows(
  df_combined,
  df_combined %>% mutate(sample = "Combined")
)

# Make sample order (8 originals + combined last)
sample_levels <- c(sort(unique(df$sample)), "Combined")
df_all$sample <- factor(df_all$sample, levels = sample_levels)

# Correlations per facet (including combined)
cor_stats <- df_all %>%
  group_by(sample) %>%
  summarise(
    rho = cor(prob, TSK1, method = "spearman", use = "complete.obs"),
    p = cor.test(prob, TSK1, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = paste0("rho = ", round(rho, 3),
                        "\nP = ", signif(p, 3)))

# Plot
p <- ggplot(df_all, aes(prob, TSK1)) +
  geom_point(alpha = 0.5, size = 1.2, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
  facet_wrap(~sample, scales = "free") +
  geom_text(
    data = cor_stats,
    aes(label = label),
    x = Inf, y = Inf,
    hjust = 1.1, vjust = 1.2,
    inherit.aes = FALSE,
    size = 3.5
  ) +
  labs(
    x = "MuT-GCNN prediction",
    y = "TSK gene signature",
    title = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

p
