library(svglite)
library(Seurat)
library(ggplot2)
library(reshape2)
library(tidyr)
library(dplyr)
library(patchwork)

#### B ####
lazyLoad('data/VISIUM/lazy_shiny_data')
# H&E
SpatialFeaturePlot(UpdateSeuratObject(Seurat_counts[[1]]),features = 'prob_scaled', ncol = 1, alpha = 0 )+ theme(legend.position = "None")
SpatialFeaturePlot(UpdateSeuratObject(Seurat_counts[[2]]),features = 'prob_scaled', ncol = 1, alpha = 0 )+ theme(legend.position = "None")
SpatialFeaturePlot(UpdateSeuratObject(Seurat_counts[[3]]),features = 'prob_scaled', ncol = 1, alpha = 0 )+ theme(legend.position = "None")
SpatialFeaturePlot(UpdateSeuratObject(Seurat_counts[[4]]),features = 'prob_scaled', ncol = 1, alpha = 0 )+ theme(legend.position = "None")



# UMAP

temp_seurat <- FindClusters(UpdateSeuratObject(Seurat_counts[[1]]), algorithm = 4, resolution = 0.15)
temp_seurat <- RunUMAP(temp_seurat, reduction = "pca", dims = 1:30)
DimPlot(temp_seurat, label = T, cols = c("#F8766D","#00BA38","#619CFF")) + theme(legend.position = "None")
SpatialDimPlot(temp_seurat,image.alpha = 0, label = T, cols = c("#F8766D","#00BA38","#619CFF")) + theme(legend.position = "None")


temp_seurat <- FindClusters(UpdateSeuratObject(Seurat_counts[[2]]), algorithm = 4, resolution = 0.10)
temp_seurat <- RunUMAP(temp_seurat, reduction = "pca", dims = 1:30)
DimPlot(temp_seurat, label = T, cols = c("#00BA38","#619CFF")) + theme(legend.position = "None")
SpatialDimPlot(temp_seurat,image.alpha = 0, label = T, cols = c("#00BA38","#619CFF")) + theme(legend.position = "None")


temp_seurat <- FindClusters(UpdateSeuratObject(Seurat_counts[[3]]), algorithm = 4, resolution = 0.1)
temp_seurat <- RunUMAP(temp_seurat, reduction = "pca", dims = 1:30)
DimPlot(temp_seurat, label = T, cols = c("#F8766D","#00BA38")) + theme(legend.position = "None")
SpatialDimPlot(temp_seurat,image.alpha = 0, label = T, cols = c("#F8766D","#00BA38")) + theme(legend.position = "None")


temp_seurat <- FindClusters(UpdateSeuratObject(Seurat_counts[[4]]), algorithm = 4, resolution = 0.15)
temp_seurat <- RunUMAP(temp_seurat, reduction = "pca", dims = 1:30)
DimPlot(temp_seurat, label = T, cols = c("#00BA38","#F8766D","#619CFF")) + theme(legend.position = "None")
SpatialDimPlot(temp_seurat,image.alpha = 0,cols = c("#00BA38","#F8766D","#619CFF"), label = T) + theme(legend.position = "None")


# SPOT-Mut prediction score
SpatialColors <- function(){
  SpColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
  return(SpColors(n = 100))
}
cols <- SpatialColors()

SpatialFeaturePlot(UpdateSeuratObject(Seurat_counts[[1]]),features = 'prob_scaled', ncol = 1, image.alpha = 0) +
  ggplot2::scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0), name = "SPOT-Mut prediction")+
  theme(legend.position = 'right',
        legend.background = element_blank())


SpatialFeaturePlot(UpdateSeuratObject(Seurat_counts[[2]]),features = 'prob_scaled', ncol = 1, image.alpha = 0) +
  ggplot2::scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0), name = "SPOT-Mut prediction")+
  theme(legend.position = 'right',
        legend.background = element_blank())


SpatialFeaturePlot(UpdateSeuratObject(Seurat_counts[[3]]),features = 'prob_scaled', ncol = 1, image.alpha = 0) +
  ggplot2::scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0), name = "SPOT-Mut prediction")+
  theme(legend.position = 'right',
        legend.background = element_blank())


SpatialFeaturePlot(UpdateSeuratObject(Seurat_counts[[4]]),features = 'prob_scaled', ncol = 1, image.alpha = 0) +
  ggplot2::scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0), name = "SPOT-Mut prediction")+
  theme(legend.position = 'right',
        legend.background = element_blank())


#### VAF distributions ####
load("data/AVITI_WES_mutations.Rdata")
## Venn Diagram
SNVs <- subset(maf_data, Variant_Type == 'SNP' & !Variant_Classification %in% c('RNA','Intron','IGR','COULD_NOT_DETERMINE',"5'Flank","3'UTR","5'UTR"))
SNVs$venn_id <- paste(SNVs$Chromosome, SNVs$Start_Position, SNVs$Tumor_Seq_Allele2, sep = '_')


metadata_table <- read.table('data/VISIUM_metadata.csv', sep = ';', header = 1)
rownames(metadata_table) <- metadata_table$NXTGNT_ID
metadata_table$Donor <- sapply(strsplit(metadata_table$Sample_Name, ' '),'[[',1)
metadata_table$pub_name <- c('TE2','TC2','TE1','TC1')
df <- data.frame(
  tumor_f = SNVs$tumor_f,
  sample = SNVs$Tumor_Sample_Barcode
)
# Replace sample names with aliases from metadata_table
df <- df %>%
  left_join(
    metadata_table %>%
      select(Sample_Name, pub_name),
    by = c("sample" = "Sample_Name")
  ) %>%
  mutate(sample = pub_name) %>%
  select(-pub_name)
# Annotation dataframe for top panel
ann_df_top <- SNVs %>%
  filter(
    Hugo_Symbol == "TP53",
    Tumor_Sample_Barcode %in% c("CRIG_13 TE 1", "CRIG_13 TC 1"),
    Variant_Classification != 'Silent'
  ) %>%
  distinct(
    tumor_f,
    cDNA_Change,
    Tumor_Sample_Barcode
  )

# Annotation dataframe for bottom panel
ann_df_bottom <- SNVs %>%
  filter(
    Hugo_Symbol == "TP53",
    Tumor_Sample_Barcode %in% c("CRIG_14 TC", "CRIG_14 TE"),
    Variant_Classification != 'Silent'
  ) %>%
  distinct(
    tumor_f,
    cDNA_Change,
    Tumor_Sample_Barcode
  )
ann_df_top$cDNA_Change <- paste("TP53",ann_df_top$cDNA_Change)
ann_df_top$Tumor_Sample_Barcode <- ifelse(ann_df_top$Tumor_Sample_Barcode == 'CRIG_13 TC 1', "TC1", "TE1")
ann_df_top <- rbind(data.frame(tumor_f = (14659/32965)/2,
                                  cDNA_Change = "TC1 cSCC expected VAF",
                                  Tumor_Sample_Barcode = 'Area4'),
                    ann_df_top)
ann_df_bottom$cDNA_Change <- paste("TP53",ann_df_bottom$cDNA_Change)
ann_df_bottom$Tumor_Sample_Barcode <- ifelse(ann_df_bottom$Tumor_Sample_Barcode == 'CRIG_14 TC', "TC2", "TE2")
ann_df_bottom <- rbind(data.frame(tumor_f = c((613/12210)/2, (3380/12210)/2, (13216/27475)/2),
                                  cDNA_Change = c("TE2 clone expected VAF","TE2 cSCC expected VAF", "TC2 cSCC expected VAF"),
                                  Tumor_Sample_Barcode = c('Area1','Area2', "Area3")),
                       ann_df_bottom)

calculated_vafs <- data.frame(calculated = c((14659/32965)/2,(13216/27475)/2,(3380/12210)/2,(613/12210)/2),
                              VAF = c(0.239,0.253,0.144,mean(0.025,0.022)))
rownames(calculated_vafs) <- c('TC1','TC2','TE1 Tumor','TE1 Clone')
# Consistent colors
sample_cols <- c(
  "TE1" = "red",
  "TC1" = "blue",
  "TE2"   = "red",
  "TC2"   = "blue",
  "Area1"   = "#28bc12",
  "Area2"   = "#12bcad",
  "Area3" = "#5353ffff",
  "Area4" = "#a512ddff"
)

# Top panel
p1 <- ggplot(
  subset(df, sample %in% c("TE1", "TC1") & tumor_f <= 0.5),
  aes(x = tumor_f, fill = sample)
) +
  geom_histogram(
    bins = 100,
    alpha = 0.5,
    position = "identity"
  ) +
  
  # TP53 marker lines
  geom_segment(
    data = ann_df_top,
    aes(
      x = tumor_f,
      xend = tumor_f,
      y = -25,
      yend = -5,
      color = Tumor_Sample_Barcode
    ),
    inherit.aes = FALSE,
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  
  # TP53 labels
  geom_text(
    data = ann_df_top,
    aes(
      x = tumor_f,
      y = -38,
      label = cDNA_Change,
      color = Tumor_Sample_Barcode
    ),
    inherit.aes = FALSE,
    size = 3.4,
    angle = 45,
    hjust = 1,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(values = sample_cols) +
  scale_color_manual(values = sample_cols) +
  
  coord_cartesian(
    xlim = c(0, 0.5),
    ylim = c(0,250),
    clip = "off"
  ) +
  
  labs(
    title = "TC1/TE1 VAF distribution",
    x = "Variant Allele Frequency",
    y = "Frequency"
  ) +
  
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

# Bottom panel
p2 <- ggplot(
  subset(df, sample %in% c("TC2", "TE2")& tumor_f <= 0.5),
  aes(x = tumor_f, fill = sample)
) +
  geom_histogram(
    bins = 100,
    alpha = 0.5,
    position = "identity"
  ) +
  geom_segment(
    data = ann_df_bottom,
    aes(
      x = tumor_f,
      xend = tumor_f,
      y = -65,
      yend = -18,
      color = Tumor_Sample_Barcode
    ),
    inherit.aes = FALSE,
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  geom_text(
    data = ann_df_bottom,
    aes(
      x = tumor_f,
      y = -100,
      label = cDNA_Change,
      color = Tumor_Sample_Barcode
    ),
    inherit.aes = FALSE,
    size = 3.4,
    angle = 45,
    hjust = 1,
    show.legend = FALSE
  ) +
  
  scale_fill_manual(values = sample_cols) +
  scale_color_manual(values = sample_cols) +
  
  coord_cartesian(
    xlim = c(0, 0.5),
    ylim = c(0, 700),
    clip = "off"
  ) +
  
  labs(
    title = "TC2/TE2 VAF distribution",
    x = "Variant Allele Frequency",
    y = "Frequency"
  ) +
  theme_bw() +
  theme(
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

# Combine plots with reduced spacing
(p1 / p2) +
  plot_layout(heights = c(1, 1)) &
  theme(plot.margin = margin(10, 10, 50, 14))

