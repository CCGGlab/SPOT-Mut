## Supplementary 2
library(Seurat)
library(RColorBrewer)
library(patchwork)

Seurat_objects <- list()
patients <- c("4","6","16","17","18","19","20")
for (Patient in patients){
  pattern <- paste0("ST_P",Patient,"_rep(\\d+)")
  Samples <- dir('data/VISIUM/',pattern = pattern)[sapply(dir('data/VISIUM/',pattern = pattern), function(x) dir.exists(paste0('data/VISIUM/',x,'/spatial')))]
  for (Sample in Samples){
    patient_folder <- paste0('data/VISIUM/',Sample)
    print(Sample)
    Seurat.query <- Load10X_Spatial(patient_folder)
    

    ## Add the GCN probabilities
    probs <- read.table(paste0(patient_folder,'/TPM_GNN_probabilities_scaled.tsv'), header = T)
    poslist <- read.csv(paste0(patient_folder,'/spatial/tissue_positions_list.csv'), header = F, row.names = 1)
    Seurat.query@images[[1]]@scale.factors[["lowres"]] <-  ScaleFactors(Seurat.query@images[[1]])$hires
    colnames(poslist)[c(2,3)]<- c('x','y')
    poslist$barcodes <- rownames(poslist)
    probs <- merge(probs, poslist, by = c('x','y'))
    rownames(probs) <- probs$barcodes
    Seurat_objects[[Sample]] <- AddMetaData(Seurat.query, metadata = probs)
  }
}

SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))
cols <- SpatialColors(n = 100)
plotlist <- list()
# svglite("results/figs/suppl_figures/S2/cSCC1.svg", width = 6, height = 6)
plotlist[[1]] <- SpatialFeaturePlot(Seurat_objects[[1]], features = 'prob', pt.size.factor = 2.5)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.title = element_blank())
# dev.off()

# svglite("results/figs/suppl_figures/S2/cSCC2.svg", width = 6, height = 6)
plotlist[[2]] <- SpatialFeaturePlot(Seurat_objects[[2]], features = 'prob', pt.size.factor =2)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/cSCC3.svg", width = 6, height = 6)
plotlist[[3]] <- SpatialFeaturePlot(Seurat_objects[[3]], features = 'prob', pt.size.factor =2)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/cSCC4.svg", width = 6, height = 6)
plotlist[[4]] <- SpatialFeaturePlot(Seurat_objects[[4]], features = 'prob', pt.size.factor =2)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/cSCC5.svg", width = 6, height = 6)
plotlist[[5]] <- SpatialFeaturePlot(Seurat_objects[[5]], features = 'prob', pt.size.factor =1.8)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/cSCC6.svg", width = 6, height = 6)
plotlist[[6]] <- SpatialFeaturePlot(Seurat_objects[[6]], features = 'prob', pt.size.factor =1.8)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/cSCC7.svg", width = 6, height = 6)
plotlist[[7]] <- SpatialFeaturePlot(Seurat_objects[[7]], features = 'prob', pt.size.factor =1.8)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/cSCC8.svg", width = 6, height = 6)
plotlist[[8]] <- SpatialFeaturePlot(Seurat_objects[[8]], features = 'prob', pt.size.factor =1.8)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/oSCC1.svg", width = 6, height = 6)
plotlist[[9]] <- SpatialFeaturePlot(Seurat_objects[[9]], features = 'prob', pt.size.factor =1.6)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/oSCC2.svg", width = 6, height = 6)
plotlist[[10]] <- SpatialFeaturePlot(Seurat_objects[[10]], features = 'prob', pt.size.factor =1.87)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/oSCC3.svg", width = 6, height = 6)
plotlist[[11]] <- SpatialFeaturePlot(Seurat_objects[[11]], features = 'prob', pt.size.factor =1.8)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/oSCC4.svg", width = 6, height = 6)
plotlist[[12]] <- SpatialFeaturePlot(Seurat_objects[[12]], features = 'prob', pt.size.factor =1.7)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/oSCC5svg", width = 6, height = 6)
plotlist[[13]] <- SpatialFeaturePlot(Seurat_objects[[13]], features = 'prob', pt.size.factor =1.8)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/oSCC6.svg", width = 6, height = 6)
plotlist[[14]] <- SpatialFeaturePlot(Seurat_objects[[14]], features = 'prob', pt.size.factor =1.8)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/oSCC7.svg", width = 6, height = 6)
plotlist[[15]] <- SpatialFeaturePlot(Seurat_objects[[15]], features = 'prob', pt.size.factor =1.7)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/oSCC8.svg", width = 6, height = 6)
plotlist[[16]] <- SpatialFeaturePlot(Seurat_objects[[16]], features = 'prob', pt.size.factor =1.8)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/oSCC9.svg", width = 6, height = 6)
plotlist[[17]] <- SpatialFeaturePlot(Seurat_objects[[17]], features = 'prob', pt.size.factor =1.7)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/oSCC10.svg", width = 6, height = 6)
plotlist[[18]] <- SpatialFeaturePlot(Seurat_objects[[18]], features = 'prob', pt.size.factor =1.8)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/oSCC11.svg", width = 6, height = 6)
plotlist[[19]] <- SpatialFeaturePlot(Seurat_objects[[19]], features = 'prob', pt.size.factor =1.9)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

# svglite("results/figs/suppl_figures/S2/oSCC12.svg", width = 6, height = 6)
plotlist[[20]] <- SpatialFeaturePlot(Seurat_objects[[20]], features = 'prob', pt.size.factor =1.9)+
  scale_fill_gradientn(limits = c(0.0,1.0), colours = cols, breaks = c(0,0.5,1.0))+
  theme(legend.position = "none")
# dev.off()

p <- wrap_plots(plotlist, nrow = 5, ncol = 4)

# ggsave(
#   filename = "Fig_S3.png",
#   plot = p,
#   width = 25,
#   height = 20,
#   units = "in"
# )

# svglite("Fig_S3.svg", width = 17, height = 17)
# p
# dev.off()