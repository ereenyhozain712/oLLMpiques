# =====================================================================
# Cannabis sativa: REM16 and FT1 expression by plant sex
# REM16 = LOC115699937, FT1 = LOC115696989
# Input: 2_Data_RNASeq_Cannabis_Sex.csv
#   - genes in rows (locus IDs in column 1), samples in columns
#   - sex encoded in the sample name suffix: _XX = female, _XY = male
# =====================================================================

# Install required packages if they are not already installed
required_packages <- c("data.table", "ggplot2")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}

# Load required packages
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ---------------------------- SETTINGS -------------------------------
# NOTE FOR THE READER: before running this script, set your working
# directory to the folder that contains the data file "2_Data_RNASeq_Cannabis_Sex.csv" as described below.
#   In RStudio: Session > Set Working Directory > To Source File Location
#   (or place this script and the CSV in the same folder).
# Outputs (figures and tables) are written to a "figures" subfolder there.
EXPR_FILE <- "2_Data_RNASeq_Cannabis_Sex.csv"
OUT_DIR   <- "figures"
Y_LABEL   <- "Normalized expression (as provided)"
GENES     <- c(REM16 = "LOC115699937", FT1 = "LOC115696989")
# ---------------------------------------------------------------------

# Stop early with a helpful message if the data file can't be found
if (!file.exists(EXPR_FILE)) {
  stop("Cannot find '", EXPR_FILE, "' in the current working directory:\n  ",
       getwd(), "\nSet the working directory to the folder containing the data ",
       "(Session > Set Working Directory > To Source File Location).")
}

dir.create(OUT_DIR, showWarnings = FALSE)
setwd_msg <- getwd(); cat("Working directory:", setwd_msg, "\n")

# ---------------------- 1. Read the data -----------------------------
dt <- fread(EXPR_FILE, data.table = FALSE, check.names = FALSE)
cat("Loaded:", nrow(dt), "genes x", ncol(dt) - 1, "samples\n")

gene_ids    <- trimws(as.character(dt[[1]]))     # column A (header cell is blank)
sample_ids  <- trimws(names(dt)[-1])

# ---------------------- 2. Structure checks --------------------------
# (a) duplicated sample names would silently break grouping
if (any(duplicated(sample_ids))) {
  warning("Duplicated sample names: ",
          paste(unique(sample_ids[duplicated(sample_ids)]), collapse = ", "))
}
# (b) duplicated gene IDs
if (any(duplicated(gene_ids))) {
  message(sum(duplicated(gene_ids)), " duplicated gene IDs in the file.")
}
# (c) all sample columns should be numeric
non_num <- names(dt)[-1][!sapply(dt[-1], is.numeric)]
if (length(non_num)) {
  message("Non-numeric columns coerced to numeric: ", paste(non_num, collapse = ", "))
}

# ---------------------- 3. Extract the two genes ---------------------
# Whole-token match, so e.g. LOC1156999370 can never match LOC115699937
find_row <- function(id) which(grepl(paste0("(^|[^0-9A-Za-z])", id, "($|[^0-9A-Za-z])"), gene_ids))

extract_gene <- function(id, nm) {
  idx <- find_row(id)
  if (length(idx) == 0) stop(nm, " (", id, ") not found in column 1.")
  if (length(idx) > 1) {
    warning(nm, " (", id, ") found in ", length(idx), " rows; averaging them.")
  }
  m <- as.matrix(dt[idx, -1, drop = FALSE])
  mode(m) <- "numeric"
  colMeans(m, na.rm = FALSE)
}

expr <- data.frame(
  sample = sample_ids,
  REM16  = unname(extract_gene(GENES[["REM16"]], "REM16")),
  FT1    = unname(extract_gene(GENES[["FT1"]],   "FT1")),
  stringsAsFactors = FALSE
)

# ---------------------- 4. Sex from sample names ---------------------
# Chromosomal sex tag at the end of the sample name; case-insensitive,
# tolerant of stray spaces. Anything else is flagged, not guessed.
tag <- toupper(sub(".*_\\s*([A-Za-z]+)\\s*$", "\\1", expr$sample))
tag[!grepl("_", expr$sample)] <- NA
expr$sex <- ifelse(tag == "XX", "Female (XX)",
                   ifelse(tag == "XY", "Male (XY)", NA))

cat("\nSample-name suffixes found:\n"); print(table(tag, useNA = "ifany"))
if (any(is.na(expr$sex))) {
  message("Unrecognized sex tag; excluded from plots: ",
          paste(expr$sample[is.na(expr$sex)], collapse = ", "))
}

# ---------------------- 5. Expression value QC -----------------------
for (nm in names(GENES)) {
  v <- expr[[nm]]
  cat(sprintf("%s (%s): %d NA, %d negative, range %.2f to %.2f\n",
              nm, GENES[[nm]], sum(is.na(v)), sum(v < 0, na.rm = TRUE),
              min(v, na.rm = TRUE), max(v, na.rm = TRUE)))
}

dat <- expr[!is.na(expr$sex), ]
dat$sex <- factor(dat$sex, levels = c("Female (XX)", "Male (XY)"))
cat("\nPlants analysed: n =", nrow(dat), "(expected 138)\n"); print(table(dat$sex))
if (nrow(dat) != 138) message("NOTE: n differs from 138 - check messages above.")

# Per-gene summary by sex, saved with the extracted data
summ <- do.call(rbind, lapply(names(GENES), function(g) {
  do.call(rbind, lapply(levels(dat$sex), function(s) {
    v <- dat[[g]][dat$sex == s]
    data.frame(gene = g, locus = GENES[[g]], sex = s, n = sum(is.finite(v)),
               median = median(v, na.rm = TRUE), mean = mean(v, na.rm = TRUE),
               sd = sd(v, na.rm = TRUE))
  }))
}))
print(summ, row.names = FALSE)
fwrite(dat[, c("sample", "sex", names(GENES))], file.path(OUT_DIR, "REM16_FT1_extracted.csv"))
fwrite(summ, file.path(OUT_DIR, "REM16_FT1_summary_by_sex.csv"))

# ---------------------- 6. Figures -----------------------------------
sex_cols <- c("Female (XX)" = "#C2185B", "Male (XY)" = "#1976D2")

theme_set(theme_bw(base_size = 13) +
            theme(panel.grid.minor = element_blank(),
                  legend.position = "none",
                  strip.background = element_rect(fill = "grey92"),
                  strip.text = element_text(face = "bold", size = 13)))

n_lab <- function(d) {
  tab <- table(d$sex)
  paste0(names(tab), "\nn = ", tab)
}

single_plot <- function(gene) {
  d <- dat[is.finite(dat[[gene]]), ]
  d$val <- d[[gene]]
  p <- tryCatch(wilcox.test(val ~ sex, data = d, exact = FALSE)$p.value,
                error = function(e) NA)
  sub <- if (is.na(p)) NULL else
    paste0("Wilcoxon rank-sum test: p = ", format.pval(p, digits = 2, eps = 1e-16))
  ggplot(d, aes(sex, val, fill = sex)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.6) +
    geom_jitter(aes(colour = sex), width = 0.15, size = 1.8, alpha = 0.7) +
    scale_fill_manual(values = sex_cols) +
    scale_colour_manual(values = sex_cols) +
    scale_x_discrete(labels = n_lab(d)) +
    labs(title = paste(gene, "expression by plant sex"),
         subtitle = sub, x = "Plant sex (chromosomal)", y = Y_LABEL)
}

p1 <- single_plot("REM16")
p2 <- single_plot("FT1")
ggsave(file.path(OUT_DIR, "Fig1_REM16_by_sex.png"), p1, width = 5, height = 5, dpi = 300)
ggsave(file.path(OUT_DIR, "Fig2_FT1_by_sex.png"),   p2, width = 5, height = 5, dpi = 300)

# Combined figure: one panel per gene, independent y-axes, sex by colour + position
long <- rbind(
  data.frame(sex = dat$sex, gene = "REM16", value = dat$REM16),
  data.frame(sex = dat$sex, gene = "FT1",   value = dat$FT1)
)
long <- long[is.finite(long$value), ]
long$gene <- factor(long$gene, levels = c("REM16", "FT1"))

p3 <- ggplot(long, aes(sex, value, fill = sex)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.6) +
  geom_jitter(aes(colour = sex), width = 0.15, size = 1.6, alpha = 0.7) +
  facet_wrap(~ gene, scales = "free_y") +
  scale_fill_manual(values = sex_cols, name = "Sex") +
  scale_colour_manual(values = sex_cols, name = "Sex") +
  scale_x_discrete(labels = n_lab(dat)) +
  labs(title = "REM16 and FT1 expression in female (XX) and male (XY) Cannabis sativa",
       subtitle = "Each point is one plant; panels use independent y-axes",
       x = "Plant sex (chromosomal)", y = Y_LABEL) +
  theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "Fig3_REM16_FT1_combined.png"), p3, width = 8.5, height = 5.5, dpi = 300)

print(p1); print(p2); print(p3)
cat("\nDone. Outputs in:", normalizePath(OUT_DIR), "\n")


# ---- Save the environment and session info ----
save.image(file.path(OUT_DIR, "REM16_FT1_workspace.RData"))   # everything in memory
saveRDS(dat, file.path(OUT_DIR, "dat_clean.rds"))             # just the cleaned data
writeLines(capture.output(sessionInfo()),
           file.path(OUT_DIR, "sessionInfo.txt"))             # R and package versions

