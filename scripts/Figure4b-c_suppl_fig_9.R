## Figure 4
## GSEA
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(svglite)
library(pathview)
library(msigdbr)


DoGSEA <- function(target_gene, analysis = 'KEGG', TERM2GENE = NA, seed = 123){
  rankdf <- read.table(paste0('data/guidedBP/',target_gene,'_gene_ranking.tsv'), header = T)
  if (target_gene == 'TP53') rankdf <- read.table('data/global_gene_ranking_values.txt')
  ranks <- rankdf[,2]
  names(ranks) <- rankdf[,1]
  
  gene_list <- sort(ranks, decreasing = TRUE)
  
  # SYMBOL -> ENTREZID mapping
  gene_map <- bitr(
    names(gene_list),
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  
  # align scores with mapped IDs
  merged <- data.frame(
    SYMBOL = names(gene_list),
    score  = as.numeric(gene_list)
  ) %>%
    inner_join(gene_map, by = "SYMBOL")
  
  ranked_vector <- merged$score
  names(ranked_vector) <- merged$ENTREZID
  ranked_vector <- sort(ranked_vector, decreasing = TRUE)
  if (is.data.frame(TERM2GENE)){
    gseadf <- GSEA(gene_list, 
                   TERM2GENE = TERM2GENE, 
                   pvalueCutoff  = 1,
                   seed = seed)
  } else {
    if (analysis == 'KEGG'){
    # KEGG GSEA
      gseadf <- gseKEGG(
        geneList     = ranked_vector,
        organism     = "hsa",
        pvalueCutoff = 1,
        pAdjustMethod = "fdr",
        verbose      = FALSE,
        seed = seed
      )
    } else if (analysis == 'GO') {
      gseadf <- gseGO(
        geneList      = ranked_vector,
        OrgDb         = org.Hs.eg.db,
        keyType       = "ENTREZID",
        ont           = "BP",       # BP, CC, MF, or ALL
        pvalueCutoff  = 1,
        pAdjustMethod = "BH",
        by            = "fgsea",
        seed = seed
      )
    } else if (analysis == 'WP'){
      gseadf <- gseWP(geneList = ranked_vector, 
                      organism = "Homo sapiens",
                      seed = seed)
      
    }
  }
  return(gseadf)
}



## Fig 4b, c        
TP53 <- DoGSEA(target_gene = "TP53", analysis = "KEGG")
dotplot(TP53, showCategory = min(sum(as.data.frame(TP53)$qvalue < 0.05), 50), x = ~ NES) + ggtitle('TP53')
gseaplot(TP53, by = "all", title = TP53$Description[1], geneSetID = 1)

## Suppl. Fig 9
KRASdf <- DoGSEA(target_gene = 'KRAS', analysis = "KEGG")
gseaplot(KRASdf, by = "all", title = KRASdf$Description[21], geneSetID = 21)
KEAP1df <- DoGSEA(target_gene = 'KEAP1', analysis = "KEGG")
gseaplot(KEAP1df, by = "all", title = KEAP1df$Description[1], geneSetID = 1)
APCdf <- DoGSEA(target_gene = 'APC', analysis = "KEGG")
gseaplot(APCdf, by = "all", title = APCdf$Description[6], geneSetID = 6)
CTNNB1df <- DoGSEA(target_gene = 'CTNNB1', analysis = "KEGG")
gseaplot(CTNNB1df, by = "all", title = CTNNB1df$Description[5], geneSetID = 5)

