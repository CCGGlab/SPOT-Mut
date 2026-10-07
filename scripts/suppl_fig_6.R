library(ggplot2)
library(ggrepel)
## Panel a ##
data <- read.delim("data/Assignment_Solution_Activities.txt", check.names = FALSE)


# Reshape data to long format
data_long <- data %>%
  pivot_longer(
    cols = -Samples,
    names_to = "Signature",
    values_to = "Count"
  ) %>%
  filter(Count > 0)
data_long$Samples <- metadata_table[data_long$Samples,"Sample_Name"]

data_long <- data_long %>%
  group_by(Samples) %>%
  mutate(TotalCount = sum(Count)) %>%
  ungroup()

data_long$Samples <- factor(
  data_long$Samples,
  levels = data_long %>%
    group_by(Samples) %>%
    summarise(TotalCount = sum(Count)) %>%
    arrange(desc(TotalCount)) %>%
    pull(Samples)
)
sigset = list('MMR_deficiency signatures' = c("6", "14", "15", "20", "21", "26", "44"),
                 'POL_deficiency signatures' = c("10a", "10b", "10c", "10d", "28"),
                 "HR_deficiency signatures" = c("3"),
                 "BER_deficiency signatures" = c("30", "36"),
                 "Chemotherapy signatures" = c("11", "25", "31", "32", "35", "86", "87", "90"),
                 'Immunosuppressants signatures' = c("32"),
                 "APOBEC_signatures" = c('2','13'),
                 "Tobacco signatures" = c("4", "29", "92"),
                 "UV signatures" = c(	"7a", "7b", "7c", "7d", "38"),
                  "Clock-like signatures"= c('1','5'),
                "Unknown" = c('33','89','19','23')
                 )
sigcols = c("#FBDD96","#9366BC","#C3AAC3","#FF7E0D","#E376C2","#BCBC21","#9ECDF7","#61646F","#FCE930","#62D59E","grey")
names(sigcols) <- names(sigset)
data_long$Signature_group <- sapply(data_long$Signature, function(x) {
  hits <- names(sigset)[sapply(sigset, function(vec) x %in% paste0('SBS',vec))]
  if (length(hits) == 0) NA else hits
})
# Plot stacked bar chart
ggplot(data_long, aes(x = Samples, y = Count, fill = Signature_group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(
    values=sigcols, name = 'Aetiology'
  )+
  labs(title = "Stacked Bar Chart of SBS Signatures in Buccal Mucosa",
       x = "Sample",
       y = "Mutation Count") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## Panel b ##
# Data
calculated_vafs <- data.frame(
  calculated = c(
    (14659 / 32965) / 2,
    (13216 / 27475) / 2,
    (3380 / 12210) / 2,
    (613 / 12210) / 2
  ),
  VAF = c(
    0.239,
    0.253,
    0.144,
    mean(c(0.025, 0.022))
  )
)

rownames(calculated_vafs) <- c(
  "TC1",
  "TC2",
  "TE1 Tumor",
  "TE1 Clone"
)

calculated_vafs$Sample <- rownames(calculated_vafs)

# Linear model (for annotation)
fit <- lm(VAF ~ calculated, data = calculated_vafs)
r2 <- summary(fit)$r.squared
eq <- sprintf(
  "y = %.3fx + %.3f\nR² = %.3f",
  coef(fit)[2],
  coef(fit)[1],
  r2
)

p <- ggplot(calculated_vafs, aes(x = calculated, y = VAF)) +
  geom_point(
    size = 4,
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 1.1
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.8,
    colour = "black",
    fill = "grey80"
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_text_repel(
    aes(label = Sample),
    size = 4,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  annotate(
    "text",
    x = min(calculated_vafs$calculated),
    y = max(calculated_vafs$VAF),
    label = eq,
    hjust = 0,
    vjust = 1,
    size = 4
  ) +
  coord_equal() +
  labs(
    x = "Calculated AF",
    y = "Observed VAF"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    panel.border = element_blank(),
    plot.margin = margin(10, 15, 10, 10)
  )


