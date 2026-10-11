# =====================================================================
# Four-method integration for one annotation level (all traits):
#   MAGMA    gene-property (gene-covar) P            [if enabled]
#   LDSC     S-LDSC (LDSC-SEG) one-sided coefficient P [if enabled for the trait]
#   scPagwas Merge_celltype_p
#   scDRS    group association MC P
# Per trait x cell type: mean -log10(P) over the available methods; bar colour =
# combination of methods with raw P < 0.05; label (#significant / #available).
# =====================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
if (exists("snakemake") && length(snakemake@log)) {
  con <- file(snakemake@log[[1]], open = "wt"); sink(con); sink(con, type = "message")
}
p <- snakemake@params
grp <- p$group; traits <- unlist(p$traits); alpha <- p$alpha
key <- function(x) gsub("[^A-Za-z0-9]", "", x)

nm <- read.csv(snakemake@input$names) %>% filter(group == grp)
cts <- nm$celltype

get_path <- function(lst, t) { v <- lst[[t]]; if (is.null(v) || !length(v) || v == "") NA else v }
rows <- lapply(traits, function(t) {
  d <- data.frame(trait = t, cell_type = cts, k = key(cts), safe = nm$safe)
  f <- get_path(p$magma, t)
  if (!is.na(f)) {
    m <- read.table(f, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
    d$MAGMA_P <- m$P[match(d$safe, m$VARIABLE)]
  } else d$MAGMA_P <- NA
  f <- get_path(p$ldsc, t)
  if (!is.na(f)) {
    l <- read.csv(f) %>% filter(group == grp)
    d$LDSC_P <- l$LDSC_P[match(d$k, key(l$celltype))]
  } else d$LDSC_P <- NA
  s <- read.csv(get_path(p$scpagwas, t))
  d$scPagwas_P <- s$pvalue[match(d$k, key(s$celltype))]
  r <- read.table(get_path(p$scdrs, t), header = TRUE, sep = "\t")
  d$scDRS_P <- r$assoc_mcp[match(d$k, key(r$group))]
  d
})
df <- bind_rows(rows) %>% select(-k, -safe)

combo <- function(P) {
  s <- c("MAGMA", "LDSC", "scPagwas", "scDRS")[!is.na(P) & P < alpha]
  if (!length(s)) "None" else if (length(s) == 4) "All four" else if (length(s) == 1) paste(s, "only") else paste(s, collapse = "+")
}
Pm <- as.matrix(df[, c("MAGMA_P", "LDSC_P", "scPagwas_P", "scDRS_P")])
df$n_sig <- rowSums(Pm < alpha, na.rm = TRUE)
df$n_avail <- rowSums(!is.na(Pm))
df$significant_methods <- apply(Pm, 1, combo)
df$mean_neglog10P <- rowMeans(-log10(pmax(Pm, 1e-300)), na.rm = TRUE)
df <- df %>% filter(n_avail > 0)
write.csv(df, snakemake@output$table, row.names = FALSE)
print(df %>% arrange(trait, desc(mean_neglog10P)) %>% mutate(across(ends_with("_P"), ~ signif(.x, 3))), row.names = FALSE)

# ---- figure (palette and layout as in the scDepBrain manuscript figures) ----
COLS <- c("All four" = "#7a0a0a",
          "MAGMA+LDSC+scPagwas" = "#cc3a36", "MAGMA+LDSC+scDRS" = "#cc3a36",
          "MAGMA+scPagwas+scDRS" = "#e85a55", "LDSC+scPagwas+scDRS" = "#e85a55",
          "MAGMA+LDSC" = "#fa8072", "MAGMA+scDRS" = "#fa8072", "LDSC+scDRS" = "#fa8072",
          "MAGMA+scPagwas" = "#ee7e60", "LDSC+scPagwas" = "#ee7e60", "scPagwas+scDRS" = "#f4a261",
          "MAGMA only" = "#6cba73", "LDSC only" = "#4c79b6", "scDRS only" = "#e88f18",
          "scPagwas only" = "#b698c5", "None" = "#b7b7b7")
ord <- df %>% group_by(cell_type) %>% summarise(m = mean(mean_neglog10P)) %>% arrange(m) %>% pull(cell_type)
df$cell_type <- factor(df$cell_type, levels = ord)
df$trait <- factor(df$trait, levels = traits)
methods_used <- c("MAGMA", "LDSC", "scPagwas", "scDRS")[colSums(!is.na(Pm)) > 0]

g <- ggplot(df, aes(mean_neglog10P, cell_type, fill = significant_methods)) +
  geom_col(width = 0.72, color = "black", linewidth = 0.25) +
  geom_text(aes(label = sprintf("%.2f", mean_neglog10P)), hjust = -0.2, size = 3) +
  geom_text(aes(x = 0, label = sprintf("(%d/%d)  ", n_sig, n_avail)), hjust = 1, size = 2.5, color = "gray40") +
  geom_vline(xintercept = -log10(alpha), linetype = 2, color = "gray45") +
  facet_wrap(~ trait, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = COLS, name = sprintf("Significant method(s)\n(raw P < %g)", alpha), drop = TRUE) +
  guides(fill = guide_legend(ncol = 3, title.position = "top")) +
  scale_x_continuous(expand = expansion(mult = c(0.18, 0.22))) +
  labs(x = expression(Mean~-log[10](P)~"across available methods"), y = NULL,
       title = sprintf("Integrated cell-type enrichment — %s", grp),
       subtitle = sprintf("%s  |  (#significant / #available methods)", paste(methods_used, collapse = " + "))) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold", size = 11), legend.position = "bottom",
        plot.title.position = "plot", legend.justification = "left",
        legend.text = element_text(size = 9), legend.title = element_text(size = 10))
w <- max(8, 2.5 + 3.6 * length(traits)); h <- 3.4 + 0.32 * length(cts)
ggsave(snakemake@output$pdf, g, width = w, height = h, limitsize = FALSE)
ggsave(snakemake@output$png, g, width = w, height = h, dpi = 200, bg = "white", limitsize = FALSE)
cat("Wrote", snakemake@output$png, "\n")
