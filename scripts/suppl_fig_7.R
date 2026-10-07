library(maftools)
library(ggplot2)
library(R.utils)
library(dplyr)
library(tidyr)
library(pheatmap)
library(VennDiagram)
library(RColorBrewer)
library(patchwork)

library(dplyr)
library(ggplot2)
library(patchwork)

load("data/AVITI_WES_mutations.Rdata")

## Filter SNVs
SNVs <- subset(
  maf_data,
  Variant_Type == "SNP" &
    !Variant_Classification %in%
    c("RNA", "Intron", "IGR", "COULD_NOT_DETERMINE",
      "5'Flank", "3'UTR", "5'UTR")
)

metadata_table <- read.table(
  "data/VISIUM/VISIUM_metadata.csv",
  sep = ";",
  header = TRUE
)

rownames(metadata_table) <- metadata_table$NXTGNT_ID
metadata_table$Donor <- sapply(
  strsplit(metadata_table$Sample_Name, " "),
  "[[",
  1
)

metadata_table$pub_name <- c("TE2", "TC2", "TE1", "TC1")


## --------------------------
## VAF dataframe
## --------------------------

df <- data.frame(
  tumor_f = SNVs$tumor_f,
  sample = SNVs$Tumor_Sample_Barcode
) %>%
  left_join(
    metadata_table %>%
      select(Sample_Name, pub_name),
    by = c("sample" = "Sample_Name")
  ) %>%
  mutate(sample = pub_name) %>%
  select(-pub_name)


## --------------------------
## TP53 annotations
## --------------------------

ann_df <- SNVs %>%
  filter(
    Hugo_Symbol == "TP53",
    Variant_Classification != "Silent",
    Tumor_Sample_Barcode %in% c(
      "CRIG_13 TC 1",
      "CRIG_13 TE 1",
      "CRIG_14 TC",
      "CRIG_14 TE"
    )
  ) %>%
  distinct(
    tumor_f,
    cDNA_Change,
    Tumor_Sample_Barcode
  ) %>%
  mutate(
    cDNA_Change = paste("TP53", cDNA_Change),
    
    sample = case_when(
      Tumor_Sample_Barcode == "CRIG_13 TC 1" ~ "TC1",
      Tumor_Sample_Barcode == "CRIG_13 TE 1" ~ "TE1",
      Tumor_Sample_Barcode == "CRIG_14 TC"   ~ "TC2",
      Tumor_Sample_Barcode == "CRIG_14 TE"   ~ "TE2"
    ),
    
    annotation_group = sample
  )


## --------------------------
## Expected VAF annotations
## --------------------------

expected_vafs <- data.frame(
  tumor_f = c(
    (14659 / 32965) / 2,   # TC1
    (13216 / 27475) / 2,   # TC2
    (3380 / 12210) / 2,    # TE2 cSCC
    (613 / 12210) / 2      # TE2 clone
  ),
  
  cDNA_Change = c(
    "expected VAF",
    "expected VAF",
    "expected VAF",
    "AK expected VAF"
  ),
  
  sample = c(
    "TC1",
    "TC2",
    "TE2",
    "TE2"
  ),
  
  annotation_group = c(
    "Area4",
    "Area3",
    "Area2",
    "Area1"
  )
)


## Combine TP53 + expected VAF annotations
ann_df <- bind_rows(
  ann_df %>%
    select(tumor_f, cDNA_Change, sample, annotation_group),
  expected_vafs
)


## --------------------------
## Colors
## --------------------------

sample_cols <- c(
  "TE1"   = "black",
  "TC1"   = "black",
  "TE2"   = "black",
  "TC2"   = "black",
  "Area1" = "#28bc12",
  "Area2" = "#12bcad",
  "Area3" = "#5353ffff",
  "Area4" = "#a512ddff"
)


## --------------------------
## Function to make one panel
## --------------------------

make_vaf_plot <- function(sample_name,
                          ymax,
                          segment_y1,
                          segment_y2,
                          text_y) {
  
  plot_df <- df %>%
    filter(
      sample == sample_name,
      tumor_f <= 0.5
    )
  
  plot_ann <- ann_df %>%
    filter(sample == sample_name)
  
  ggplot(
    plot_df,
    aes(x = tumor_f, fill = sample)
  ) +
    
    geom_histogram(
      bins = 100,
      # alpha = 0.5,
      position = "identity"
    ) +
    
    ## Annotation lines
    geom_segment(
      data = plot_ann,
      aes(
        x = tumor_f,
        xend = tumor_f,
        y = segment_y1,
        yend = segment_y2,
        color = annotation_group
      ),
      inherit.aes = FALSE,
      linewidth = 0.7,
      show.legend = FALSE
    ) +
    
    ## Annotation labels
    geom_text(
      data = plot_ann,
      aes(
        x = tumor_f,
        y = text_y,
        label = cDNA_Change,
        color = annotation_group
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
      ylim = c(0, ymax),
      clip = "off"
    ) +
    
    labs(
      # title = sample_name,
      x = "Variant Allele Frequency",
      y = "SNV count"
    ) +
    
    theme_bw() +
    
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 55, 14)
    )
}


## --------------------------
## Generate four panels
## --------------------------

p_TC1 <- make_vaf_plot(
  sample_name = "TC1",
  ymax = 250,
  segment_y1 = -15,
  segment_y2 = -7,
  text_y = -15
)

p_TE1 <- make_vaf_plot(
  sample_name = "TE1",
  ymax = 250,
  segment_y1 = -15,
  segment_y2 = -7,
  text_y = -15
)

p_TC2 <- make_vaf_plot(
  sample_name = "TC2",
  ymax = 700,
  segment_y1 = -45,
  segment_y2 = -20,
  text_y = -50
)

p_TE2 <- make_vaf_plot(
  sample_name = "TE2",
  ymax = 700,
  segment_y1 = -45,
  segment_y2 = -20,
  text_y = -50
)

p_TC1 / p_TE1 / p_TC2 / p_TE2
## --------------------------
## Combine into 4 panels
## --------------------------

(p_TC1 | p_TE1) /
  (p_TC2 | p_TE2)

ggsave(plot = p_TC1, filename = 'results/figs/Fig 6/TC1_VAF.pdf', width = 15*0.28, height = 9*0.4)
ggsave(plot = p_TE1, filename = 'results/figs/Fig 6/TE1_VAF.pdf', width = 15*0.28, height = 9*0.4)
ggsave(plot = p_TC2, filename = 'results/figs/Fig 6/TC2_VAF.pdf', width = 15*0.28, height = 9*0.4)
ggsave(plot = p_TE2, filename = 'results/figs/Fig 6/TE2_VAF.pdf', width = 15*0.28, height = 9*0.4)

plot_stacked_mutations <- function(counts, altcol = "#E45756", altname = 'Alternative'){
  
  rownames(counts) <- gsub("_", " ", rownames(counts))
  
  plot_data <- data.frame(
    sample = factor(rep(rownames(counts), each = 2), levels = rownames(counts)),
    allele = factor(
      rep(c("ref", "alt"), times = nrow(counts)),
      levels = c("ref", "alt")
    ),
    count = as.vector(t(counts))
  )
  
  labels <- data.frame(
    sample = rownames(counts),
    alt_frac = paste0(
      counts[, "alt"], "/", counts[, "ref"] + counts[, "alt"]
    ),
    VAF = round(
      counts[, "alt"] / (counts[, "alt"] + counts[, "ref"]),
      digits = 2
    ),
    y = counts[, "alt"] / (counts[, "alt"] + counts[, "ref"])
  )
  
  # Add VAF to alt rows
  plot_data$VAF <- rep(
    counts[, "alt"] / (counts[, "alt"] + counts[, "ref"]),
    each = 2
  )
  
  return(
    ggplot(
      subset(plot_data, allele == "alt"),
      aes(x = sample, y = VAF, fill = allele)
    ) +
      geom_col(
        width = 0.55,
        colour = altcol
      ) +
      geom_text(
        data = labels,
        aes(x = sample, y = VAF, label = alt_frac),
        inherit.aes = FALSE,
        vjust = -0.5,
        fontface = "bold"
      ) +
      scale_y_continuous(
        expand = expansion(mult = c(0, 0.1))
      ) +
      scale_fill_manual(
        values = c(alt = altcol),
        labels = c(alt = altname)
      )+
      labs(
        x = "Sample",
        y = "Mutated fraction",
        fill = "Mutation"
      ) +
      coord_cartesian(ylim = c(0, 0.5)) +
      theme_bw()+
      theme(
        panel.grid.minor = element_blank()
      )
  )
}
  

TC1_splice <- rbind(
  bulk = c(ref = 240, alt = 68),
  tumor  = c(ref = 256, alt = 77),
  epi  = c(ref = 153, alt = 20),
  normal = c(ref = 252, alt = 2)
)
TC1_read_count <- plot_stacked_mutations(TC1_splice, altname = "c.672G>A", altcol = "#30b31d")
ggsave(plot = TC1_read_count, filename = 'results/figs/Fig 6/TC1_read_count.pdf', width = 15*0.28, height = 9*0.4)


## TC2 chr17:7670716; c.e10+1G>T
# Combine samples into a matrix
TC2_splice <- rbind(
  bulk = c(ref = 98, alt = 31),
  tumor_1 = c(ref = 27, alt = 9),
  tumor_2 = c(ref = 33, alt = 10),
  normal = c(ref = 27, alt = 2)
)
TC2_read_count <- plot_stacked_mutations(TC2_splice, altcol = "#0072B2", altname = "c.e10+1G>T")
ggsave(plot = TC2_read_count, filename = 'results/figs/Fig 6/TC2_read_count.pdf', width = 15*0.28, height = 9*0.4)

## TE2 chr17:7670716; c.e10+1G>T
# Combine samples into a matrix

TE2_splice <- data.frame(
  sample = c("bulk", "tumor", "epi", "normal"),
  ref = c(258, 36, 39, 11),
  alt = c(40, 13, 3, 0),
  mutation = "c.e10+1G>T"
)

# plot_stacked_mutations(TE2_splice) + ggtitle("TE2 c.e10+1G>T")

## TE2 chr17:7673764; c.856G>A
# Combine samples into a matrix
TE2_mis1 <- data.frame(
  sample = c("bulk", "tumor", "epi", "normal"),
  ref = c(376, 67, 45, 8),
  alt = c(9, 0, 4, 0),
  mutation = "c.856G>A"
)

## TE2 chr17:7673776; c.844C>T
# Combine samples into a matrix
TE2_mis2 <- data.frame(
  sample = c("bulk", "tumor", "epi", "normal"),
  ref = c(359, 56, 46, 8),
  alt = c(7, 0, 4, 0),
  mutation = "c.844C>T"
)
TE2matrix <- bind_rows(
  TE2_splice,
  TE2_mis1,
  TE2_mis2
)

# Calculate VAF
TE2_plot_data <- TE2matrix %>%
  mutate(
    frequency = alt / (ref + alt),
    sample = factor(
      sample,
      levels = c("bulk", "tumor", "epi", "normal")
    )
  )
TE2_plot_data <- TE2matrix %>%
  mutate(
    frequency = alt / (ref + alt),
    label = paste0(alt, "/", alt + ref),
    sample = factor(
      sample,
      levels = c("bulk", "tumor", "epi", "normal")
    )
  )

mutation_cols <- c(
  "c.e10+1G>T" = "#0072B2",
  "c.856G>A"   = "#D55E00",
  "c.844C>T"   = "#009E73"
)
TE2_read_count <- ggplot(
  TE2_plot_data,
  aes(x = sample, y = frequency, fill = mutation)
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = label, group = mutation),
    position = position_dodge(width = 0.8),
    vjust = -0.4,
    size = 3.5
  ) +
  scale_fill_manual(values = mutation_cols) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    x = "Sample",
    y = "Mutated fraction",
    fill = "Mutation"
  ) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 0.5)) +
  theme(
    panel.grid.minor = element_blank()
  )
ggsave(plot = TE2_read_count, filename = 'results/figs/Fig 6/TE2_read_count.pdf', width = 15*0.28, height = 9*0.4)


