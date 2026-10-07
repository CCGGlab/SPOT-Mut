library(ggplot2)
library(ggrepel)
library(viridis)
library(scales)
library(tidyverse)
library(yardstick)
library(patchwork)

metrics <- read.csv("data/GCN_test_metrics.csv")
colnames(metrics) <- gsub('test_','',colnames(metrics))

curve_data <- read.csv("data/GCN_test_curve_data.csv")

curve_data$label <- as.integer(as.character(curve_data$label))

# Add genes here that you want to label manually
# manual_highlight_genes <- character(0)
manual_highlight_genes <- c("TP53", "KEAP1", "KRAS", "BRAF", "CTNNB1","APC")
highlight_genes <- unique(manual_highlight_genes)

highlight_colours <- setNames(
  grDevices::colorRampPalette(c(
    "#5E4FA2",
    "#3288BD",
    "#1B9E77",
    "#D9A400",
    "#E66101",
    "#C51B7D"
  ))(length(highlight_genes)),
  highlight_genes
)

manual_highlight_genes <- intersect(
  manual_highlight_genes,
  metrics$gene
)

highlight_genes <- unique(c(manual_highlight_genes))
highlight_data <- metrics[metrics$gene %in% highlight_genes,]

genes <- unique(curve_data$gene)

roc_list <- vector("list", length(genes))
prc_list <- vector("list", length(genes))

auc_data <- data.frame(
  gene = genes,
  AUROC = NA_real_,
  AUPRC = NA_real_
)

for (i in seq_along(genes)) {
  gene_name <- genes[i]
  
  gene_data <- curve_data[curve_data$gene == gene_name, ]
  
  gene_input <- data.frame(
    label = factor(
      gene_data$label,
      levels = c(0, 1)
    ),
    prediction_score = gene_data$prediction_score
  )
  
  gene_roc <- yardstick::roc_curve(
    gene_input,
    label,
    prediction_score,
    event_level = "second"
  )
  
  gene_prc <- yardstick::pr_curve(
    gene_input,
    label,
    prediction_score,
    event_level = "second"
  )
  
  gene_roc$gene <- gene_name
  gene_prc$gene <- gene_name
  
  roc_list[[i]] <- gene_roc
  prc_list[[i]] <- gene_prc
  
  auc_data$AUROC[i] <- yardstick::roc_auc_vec(
    gene_input$label,
    gene_input$prediction_score,
    event_level = "second"
  )
  
  auc_data$AUPRC[i] <- yardstick::pr_auc_vec(
    gene_input$label,
    gene_input$prediction_score,
    event_level = "second"
  )
}

roc_data <- do.call(rbind, roc_list)
prc_data <- do.call(rbind, prc_list)

roc_data$false_positive_rate <- 1 - roc_data$specificity

roc_data$AUROC <- auc_data$AUROC[
  match(roc_data$gene, auc_data$gene)
]

prc_data$AUPRC <- auc_data$AUPRC[
  match(prc_data$gene, auc_data$gene)
]

prc_data <- prc_data[
  is.finite(prc_data$recall) &
    is.finite(prc_data$precision),
]


# =============================================================================
# Select curve positions for gene labels
# =============================================================================

find_curve_label_points <- function(
    data,
    selected_genes,
    x_column,
    target_x
) {
  if (length(selected_genes) == 0) {
    return(data[0, , drop = FALSE])
  }
  
  label_list <- vector(
    "list",
    length(selected_genes)
  )
  
  for (i in seq_along(selected_genes)) {
    gene_name <- selected_genes[i]
    
    gene_curve <- data[
      data$gene == gene_name,
      ,
      drop = FALSE
    ]
    
    nearest_point <- which.min(
      abs(gene_curve[[x_column]] - target_x)
    )
    
    label_list[[i]] <- gene_curve[
      nearest_point,
      ,
      drop = FALSE
    ]
  }
  
  label_points <- do.call(rbind, label_list)
  rownames(label_points) <- NULL
  
  return(label_points)
}

roc_label_data <- find_curve_label_points(
  roc_data,
  manual_highlight_genes,
  "false_positive_rate",
  0.05
)

prc_label_data <- find_curve_label_points(
  prc_data,
  manual_highlight_genes,
  "recall",
  0.95
)


# =============================================================================
# ROC curves
# =============================================================================

roc_plot <- ggplot(
  roc_data,
  aes(
    x = false_positive_rate,
    y = sensitivity,
    group = gene
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    colour = "grey55",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  # All genes in light grey
  geom_path(
    colour = "grey80",
    linewidth = 0.65,
    alpha = 0.75
  ) +
  
  # Highlight genes using the original AUROC colour scale
  geom_path(
    data = roc_data[roc_data$gene %in% highlight_genes, ],
    aes(colour = gene),
    linewidth = 0.9,
    alpha = 1
  ) +
  
  geom_text_repel(
    data = roc_label_data,
    aes(
      x = false_positive_rate,
      y = sensitivity,
      label = gene
    ),
    inherit.aes = FALSE,
    colour = "black",
    seed = 1001,
    size = 2.8,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  scale_colour_manual(
    values = highlight_colours,
    breaks = highlight_genes,
    name = "Gene"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  coord_equal(
    clip = "off"
  ) +
  labs(subtitle = paste0("Median AUROC = ",sprintf("%.2f",median(auc_data$AUROC, na.rm = TRUE))),
    x = "False-positive rate",
    y = "True-positive rate"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(5.5, 15, 5.5, 5.5)
  )


# =============================================================================
# Precision-recall curves
# =============================================================================

prc_plot <- ggplot(
  prc_data,
  aes(
    x = recall,
    y = precision,
    group = gene
  )
) +
  geom_path(
    linewidth = 0.65,
    alpha = 0.75,
    colour = "grey80"
  ) +
  geom_path(
    data = prc_data[prc_data$gene %in% highlight_genes, ],
    aes(colour = gene),
    linewidth = 0.9,
    alpha = 1
  ) +
  geom_text_repel(
    data = prc_label_data,
    aes(
      x = recall,
      y = precision,
      label = gene
    ),
    inherit.aes = FALSE,
    colour = "black",
    seed = 1001,
    size = 2.8,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  scale_colour_manual(
    values = highlight_colours,
    breaks = highlight_genes,
    name = "Gene"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  coord_equal(
    clip = "off"
  ) +
  labs(subtitle = paste0(
      "Median AUPRC = ",
      sprintf(
        "%.2f",
        median(auc_data$AUPRC, na.rm = TRUE)
      )
    ),
    x = "Recall",
    y = "Precision"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(5.5, 15, 5.5, 5.5)
  )

# =============================================================================
# CCLE ROC curves
# =============================================================================

roc_data_CCLE <- read_csv("data/CCLE_ROC_curves.csv")
auc_data_CCLE <- read_tsv("data/CCLE_gene_metrics.tsv")

roc_highlight_genes <- c('BRAF','TP53','KEAP1', 'KRAS', 'CTNNB1', 'APC')

roc_label_x <- 0.1
roc_label_data_CCLE <- roc_data_CCLE %>%
  filter(gene %in% roc_highlight_genes) %>%
  group_by(gene) %>%
  slice_min(
    order_by = abs(false_positive_rate - roc_label_x),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

roc_plot_CCLE <- ggplot(
  roc_data_CCLE,
  aes(
    x = false_positive_rate,
    y = sensitivity,
    group = gene
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    colour = "grey55",
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_path(
    linewidth = 0.65,
    alpha = 0.75,
    colour = 'grey80'
  ) +
  geom_path(
    data = roc_data_CCLE[roc_data_CCLE$gene %in% manual_highlight_genes, ],
    aes(colour = gene),
    linewidth = 0.9,
    alpha = 1
  ) +
  geom_text_repel(
    data = roc_label_data_CCLE,
    aes(
      x = false_positive_rate,
      y = sensitivity,
      label = gene
    ),
    inherit.aes = FALSE,
    colour = "black",
    seed = 1001,
    size = 2.8,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  scale_colour_manual(
    values = highlight_colours,
    breaks = highlight_genes,
    name = "Gene"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  coord_equal(
    clip = "off"
  ) +
  labs(
    subtitle = paste0(
      "Median AUROC = ",
      sprintf("%.2f", median(auc_data_CCLE$AUROC, na.rm = TRUE))
    ),
    x = "False-positive rate",
    y = "True-positive rate"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(5.5, 15, 5.5, 5.5)
  )


# =============================================================================
# CCLE Precision-recall curves
# =============================================================================
pr_data_CCLE  <- read_csv("/home/tom/projects/MutGNN/results/CGC/noVal_noWeight/CCLE_PR_curves.csv")

pr_highlight_genes <- roc_highlight_genes
# Target position on x-axis for labels
pr_label_x <- 0.9

pr_label_data <- pr_data_CCLE %>%
  filter(gene %in% pr_highlight_genes) %>%
  group_by(gene) %>%
  slice_min(
    order_by = abs(recall - pr_label_x),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

pr_plot <- ggplot(
  pr_data_CCLE,
  aes(
    x = recall,
    y = precision,
    group = gene
  )
) +
  geom_path(
    linewidth = 0.65,
    alpha = 0.75,
    colour = 'grey80'
  ) +
  geom_path(
    data = pr_data_CCLE[pr_data_CCLE$gene %in% manual_highlight_genes, ],
    aes(colour = AUPRC),
    linewidth = 0.9,
    alpha = 1
  ) +
  geom_text_repel(
    data = pr_label_data,
    aes(
      x = recall,
      y = precision,
      label = gene
    ),
    inherit.aes = FALSE,
    colour = "black",
    seed = 1001,
    size = 2.8,
    box.padding = 0.3,
    point.padding = 0.2,
    min.segment.length = 0,
    max.overlaps = Inf
  ) +
  scale_colour_viridis_c(
    option = "C",
    limits = c(0, 1),
    oob = scales::squish,
    name = "AUPRC"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    expand = c(0, 0)
  ) +
  coord_equal(
    clip = "off"
  ) +
  labs(
    subtitle = paste0(
      "Median AUPRC = ",
      sprintf("%.2f", median(auc_data_CCLE$AUPRC, na.rm = TRUE))
    ),
    x = "Recall",
    y = "Precision"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    plot.margin = margin(5.5, 15, 5.5, 5.5)
  )

# combine
# ---------------------

layout <- "
AB
CD
"
combined_plot <-
  roc_plot +
  prc_plot +
  roc_plot_CCLE +
  pr_plot +
  plot_layout(
    design = layout,
    widths = c(1, 1)
  )

combined_plot
