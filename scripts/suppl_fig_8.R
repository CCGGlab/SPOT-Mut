library(ggplot2)
library(data.table)

# Load required data
load("data/gene_cancer_table.RData")

CGC <- read.csv("data/cancer_gene_census.csv")

# -------------------------------------------------------------------------
# Gene ordering
# -------------------------------------------------------------------------

genes <- intersect(
  names(result_F1),
  names(result_BA)
)

ba_data <- data.frame(
  gene = genes,
  model_BA = vapply(
    result_BA[genes],
    function(x) x$model_f1,
    numeric(1)
  )
)

gene_levels <- ba_data$gene[
  order(ba_data$model_BA)
]

shared_y_limits <- c(
  0.4,
  length(gene_levels) + 0.6
)

gene_guidelines <- seq(
  0.5,
  length(gene_levels) + 0.5,
  by = 1
)

guide_color <- "grey80"
guide_width <- 0.3

tier_colours <- c(
  "1" = "#9E0142",
  "2" = "#5E4FA2",
  "Not listed" = "#707070"
)

# -------------------------------------------------------------------------
# Prepare inclusion data
# -------------------------------------------------------------------------

sample_data <- rbindlist(
  gene_cancer_tables,
  idcol = "gene"
)

# Summary per gene
gene_summary <- sample_data[, .(
  total_size = sum(sample_size)
), by = gene]

gene_summary <- gene_summary[
  match(gene_levels, gene_summary$gene),
]

gene_label_data <- copy(gene_summary)

tier_index <- match(
  gene_label_data$gene,
  CGC$Gene.Symbol
)

gene_label_data$tier <- as.character(
  CGC$Tier[tier_index]
)

gene_label_data$tier[
  is.na(gene_label_data$tier)
] <- "Not listed"

gene_label_data$y_position <- match(
  gene_label_data$gene,
  gene_levels
)

gene_label_data$label <- paste0(
  gene_label_data$gene,
  " (n = ",
  gene_label_data$total_size,
  ")"
)

# Summary per cancer type
cancer_summary <- sample_data[
  ,
  .(
    n_inclusions = uniqueN(gene),
    sample_size = sample_size[1]
  ),
  by = cancer_type
]

cancer_summary <- cancer_summary[
  order(
    n_inclusions,
    sample_size,
    decreasing = TRUE
  ),
]

cancer_levels <- cancer_summary$cancer_type

cancer_labels <- paste0(
  cancer_summary$cancer_type,
  " (n = ",
  cancer_summary$sample_size,
  ")"
)

# -------------------------------------------------------------------------
# Complete gene × cancer-type grid
# -------------------------------------------------------------------------

sample_data[, mutation_proportion := mutated / sample_size]

grid_data <- expand.grid(
  gene = gene_levels,
  cancer_type = cancer_levels,
  stringsAsFactors = FALSE
)

combination_key <- paste(
  sample_data$gene,
  sample_data$cancer_type
)

grid_key <- paste(
  grid_data$gene,
  grid_data$cancer_type
)

grid_data$mutation_proportion <- sample_data$mutation_proportion[
  match(grid_key, combination_key)
]

grid_data$proportion_label <- ifelse(
  is.na(grid_data$mutation_proportion),
  "",
  format(
    round(grid_data$mutation_proportion, 2),
    nsmall = 1,
    trim = TRUE
  )
)

grid_data$x_position <- match(
  grid_data$cancer_type,
  cancer_levels
)

grid_data$y_position <- match(
  grid_data$gene,
  gene_levels
)

# -------------------------------------------------------------------------
# Plot
# -------------------------------------------------------------------------

p_inclusions <- ggplot() +

  geom_hline(
    yintercept = gene_guidelines,
    colour = guide_color,
    linewidth = guide_width
  ) +

  geom_tile(
    data = grid_data,
    aes(
      x = x_position,
      y = y_position,
      fill = mutation_proportion
    ),
    colour = "white",
    linewidth = 0.4,
    width = 0.9,
    height = 0.9
  ) +

  geom_text(
    data = grid_data[
      !is.na(grid_data$mutation_proportion) &
        grid_data$mutation_proportion < 0.5,
    ],
    aes(
      x = x_position,
      y = y_position,
      label = proportion_label
    ),
    colour = "white",
    size = 2.5
  ) +

  geom_text(
    data = grid_data[
      !is.na(grid_data$mutation_proportion) &
        grid_data$mutation_proportion >= 0.5,
    ],
    aes(
      x = x_position,
      y = y_position,
      label = proportion_label
    ),
    colour = "black",
    size = 2.5
  ) +

  coord_cartesian(
    xlim = c(0.5, length(cancer_levels) + 0.5),
    clip = "off"
  ) +

  geom_text(
    data = gene_label_data,
    aes(
      x = 0.35,
      y = y_position,
      label = label,
      colour = tier
    ),
    hjust = 1,
    fontface = "italic",
    size = 3.5
  ) +

  scale_fill_viridis_c(
    name = "Mutated samples",
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    option = "C",
    na.value = "white"
  ) +

  scale_colour_manual(
    name = "CGC tier",
    values = tier_colours,
    na.value = "#707070"
  ) +

  scale_x_continuous(
    name = "Cancer type",
    breaks = seq_along(cancer_levels),
    labels = cancer_labels,
    expand = expansion(add = 0.5)
  ) +

  scale_y_continuous(
    name = NULL,
    breaks = NULL,
    limits = shared_y_limits,
    expand = c(0, 0)
  ) +

  guides(
    colour = guide_legend(order = 1),
    fill = guide_colourbar(order = 2)
  ) +

  theme_minimal() +

  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "top",
    plot.margin = margin(10, 10, 10, 100)
  )

p_inclusions