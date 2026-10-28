# =====================================================================
# Cannabis sativa: REM16 and FT1 expression by plant sex
# REM16 = LOC115699937, FT1 = LOC115696989
#
# Run from any directory:
#   Rscript /path/to/oLLMpiques/code/sex_gene_expression.R
#
# Outputs:
#   data/REM16_FT1_extracted.csv
#   data/REM16_FT1_summary_by_sex.csv
#   data/sessionInfo.txt
#   figures/Fig1_REM16_by_sex.png
#   figures/Fig2_FT1_by_sex.png
#   figures/Fig3_REM16_FT1_combined.png
# =====================================================================

# ----------------------- 1. Locate project ----------------------------

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)

if (length(file_arg) > 0) {
  script_path <- normalizePath(
    sub("^--file=", "", file_arg[1]),
    mustWork = TRUE
  )
  script_dir <- dirname(script_path)

  project_dir <- if (basename(script_dir) == "code") {
    dirname(script_dir)
  } else {
    script_dir
  }
} else {
  project_dir <- getwd()
}

project_dir <- normalizePath(project_dir, mustWork = TRUE)

EXPR_FILE <- file.path(
  project_dir, "data", "2_Data_RNASeq_Cannabis_Sex.csv"
)

DATA_DIR <- file.path(project_dir, "data")
FIGURES_DIR <- file.path(project_dir, "figures")

Y_LABEL <- "Normalized expression (as provided)"
GENES <- c(REM16 = "LOC115699937", FT1 = "LOC115696989")

dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Project directory:", project_dir, "\n")
cat("Input file:", EXPR_FILE, "\n")

if (!file.exists(EXPR_FILE)) {
  stop(
    "Cannot find input file:\n", EXPR_FILE,
    "\nCheck that the CSV is in the repository's data/ folder."
  )
}

# ----------------------- 2. Install packages --------------------------

required_packages <- c("data.table", "ggplot2")

# Use a writable user-specific R library.
user_library <- path.expand("~/R/library")
dir.create(user_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(user_library, .libPaths())))

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  cat(
    "Installing missing packages:",
    paste(missing_packages, collapse = ", "),
    "\n"
  )

  install.packages(
    missing_packages,
    lib = user_library,
    repos = "https://cloud.r-project.org"
  )
}

still_missing <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(still_missing) > 0) {
  stop(
    "Could not install/load these packages: ",
    paste(still_missing, collapse = ", "),
    "\nOn Ubuntu, if compilation tools are missing, run:\n",
    "sudo apt update && sudo apt install build-essential"
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ----------------------- 3. Read data --------------------------------

dt <- fread(
  EXPR_FILE,
  data.table = FALSE,
  check.names = FALSE
)

if (ncol(dt) < 2) {
  stop("The input file must contain a gene-ID column and sample columns.")
}

cat(
  "Loaded:", nrow(dt), "genes x",
  ncol(dt) - 1, "samples\n"
)

gene_ids <- trimws(as.character(dt[[1]]))
sample_ids <- trimws(names(dt)[-1])

if (anyNA(gene_ids) || any(!nzchar(gene_ids))) {
  warning("Some gene IDs in column 1 are empty or missing.")
}

if (anyNA(sample_ids) || any(!nzchar(sample_ids))) {
  stop("One or more sample columns have missing/empty column names.")
}

if (anyDuplicated(sample_ids)) {
  warning(
    "Duplicated sample names: ",
    paste(unique(sample_ids[duplicated(sample_ids)]), collapse = ", ")
  )
}

if (anyDuplicated(gene_ids)) {
  message(
    sum(duplicated(gene_ids)),
    " duplicated gene IDs in the file."
  )
}

# Coerce sample columns to numeric, as in the original script.
non_num <- names(dt)[-1][
  !vapply(dt[-1], is.numeric, logical(1))
]

if (length(non_num) > 0) {
  message(
    "Non-numeric columns coerced to numeric: ",
    paste(non_num, collapse = ", ")
  )

  dt[-1] <- lapply(
    dt[-1],
    function(x) suppressWarnings(as.numeric(x))
  )
}

# ----------------------- 4. Extract target genes ---------------------

# Match each exact locus ID as a whole token.
# Other gene names do not need to begin with "LOC".
find_row <- function(id) {
  which(
    grepl(
      paste0("(^|[^0-9A-Za-z])", id, "($|[^0-9A-Za-z])"),
      gene_ids
    )
  )
}

extract_gene <- function(id, gene_name) {
  idx <- find_row(id)

  if (length(idx) == 0) {
    stop(
      gene_name, " (", id, ") not found in the first column.\n",
      "Check the exact gene/locus ID in the input CSV."
    )
  }

  if (length(idx) > 1) {
    warning(
      gene_name, " (", id, ") found in ", length(idx),
      " rows; averaging them."
    )
  }

  m <- as.matrix(dt[idx, -1, drop = FALSE])
  mode(m) <- "numeric"

  colMeans(m, na.rm = FALSE)
}

expr <- data.frame(
  sample = sample_ids,
  REM16 = unname(extract_gene(GENES[["REM16"]], "REM16")),
  FT1 = unname(extract_gene(GENES[["FT1"]], "FT1")),
  stringsAsFactors = FALSE
)

# ----------------------- 5. Determine plant sex ----------------------

tag <- toupper(
  sub(".*_\\s*([A-Za-z]+)\\s*$", "\\1", expr$sample)
)

tag[!grepl("_", expr$sample)] <- NA_character_

expr$sex <- ifelse(
  tag == "XX", "Female (XX)",
  ifelse(tag == "XY", "Male (XY)", NA_character_)
)

cat("\nSample-name suffixes found:\n")
print(table(tag, useNA = "ifany"))

if (anyNA(expr$sex)) {
  message(
    "Unrecognized sex tag; excluded from plots: ",
    paste(expr$sample[is.na(expr$sex)], collapse = ", ")
  )
}

# ----------------------- 6. Expression QC ----------------------------

for (gene_name in names(GENES)) {
  v <- expr[[gene_name]]
  finite_v <- v[is.finite(v)]

  cat(sprintf(
    "%s (%s): %d NA, %d negative, range %s to %s\n",
    gene_name,
    GENES[[gene_name]],
    sum(is.na(v)),
    sum(v < 0, na.rm = TRUE),
    if (length(finite_v)) {
      sprintf("%.2f", min(finite_v))
    } else {
      "NA"
    },
    if (length(finite_v)) {
      sprintf("%.2f", max(finite_v))
    } else {
      "NA"
    }
  ))
}

dat <- expr[!is.na(expr$sex), , drop = FALSE]

dat$sex <- factor(
  dat$sex,
  levels = c("Female (XX)", "Male (XY)")
)

cat("\nPlants analysed: n =", nrow(dat), "(expected 138)\n")
print(table(dat$sex))

if (nrow(dat) != 138) {
  message(
    "NOTE: n differs from 138 - check the sample names and input data."
  )
}

# ----------------------- 7. Summary tables ---------------------------

summ <- do.call(
  rbind,
  lapply(names(GENES), function(gene_name) {
    do.call(
      rbind,
      lapply(levels(dat$sex), function(sex_name) {
        v <- dat[[gene_name]][dat$sex == sex_name]

        data.frame(
          gene = gene_name,
          locus = GENES[[gene_name]],
          sex = sex_name,
          n = sum(is.finite(v)),
          median = median(v, na.rm = TRUE),
          mean = mean(v, na.rm = TRUE),
          sd = sd(v, na.rm = TRUE)
        )
      })
    )
  })
)

print(summ, row.names = FALSE)

# Save tabular results in data/, not figures/.
fwrite(
  dat[, c("sample", "sex", names(GENES))],
  file.path(DATA_DIR, "REM16_FT1_extracted.csv")
)

fwrite(
  summ,
  file.path(DATA_DIR, "REM16_FT1_summary_by_sex.csv")
)

# ----------------------- 8. Figures ----------------------------------

sex_cols <- c(
  "Female (XX)" = "#C2185B",
  "Male (XY)" = "#1976D2"
)

theme_set(
  theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "none",
      strip.background = element_rect(fill = "grey92"),
      strip.text = element_text(face = "bold", size = 13)
    )
)

n_lab <- function(d) {
  tab <- table(d$sex)
  paste0(names(tab), "\nn = ", as.integer(tab))
}

single_plot <- function(gene_name) {
  d <- dat[is.finite(dat[[gene_name]]), , drop = FALSE]
  d$val <- d[[gene_name]]

  p_value <- tryCatch(
    wilcox.test(
      val ~ sex,
      data = d,
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )

  subtitle <- if (is.na(p_value)) {
    NULL
  } else {
    paste0(
      "Wilcoxon rank-sum test: p = ",
      format.pval(p_value, digits = 2, eps = 1e-16)
    )
  }

  ggplot(d, aes(sex, val, fill = sex)) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.55,
      alpha = 0.6
    ) +
    geom_jitter(
      aes(colour = sex),
      width = 0.15,
      size = 1.8,
      alpha = 0.7
    ) +
    scale_fill_manual(values = sex_cols) +
    scale_colour_manual(values = sex_cols) +
    scale_x_discrete(labels = n_lab(d)) +
    labs(
      title = paste(gene_name, "expression by plant sex"),
      subtitle = subtitle,
      x = "Plant sex (chromosomal)",
      y = Y_LABEL
    )
}

p1 <- single_plot("REM16")
p2 <- single_plot("FT1")

# Save each plot directly to PNG. Do not print plots to the default
# graphics device, which can create Rplots.pdf in non-interactive R.
ggsave(
  filename = file.path(FIGURES_DIR, "Fig1_REM16_by_sex.png"),
  plot = p1,
  width = 5,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path(FIGURES_DIR, "Fig2_FT1_by_sex.png"),
  plot = p2,
  width = 5,
  height = 5,
  dpi = 300
)

long <- rbind(
  data.frame(
    sex = dat$sex,
    gene = "REM16",
    value = dat$REM16
  ),
  data.frame(
    sex = dat$sex,
    gene = "FT1",
    value = dat$FT1
  )
)

long <- long[is.finite(long$value), , drop = FALSE]
long$gene <- factor(long$gene, levels = c("REM16", "FT1"))

p3 <- ggplot(long, aes(sex, value, fill = sex)) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.55,
    alpha = 0.6
  ) +
  geom_jitter(
    aes(colour = sex),
    width = 0.15,
    size = 1.6,
    alpha = 0.7
  ) +
  facet_wrap(~gene, scales = "free_y") +
  scale_fill_manual(values = sex_cols, name = "Sex") +
  scale_colour_manual(values = sex_cols, name = "Sex") +
  scale_x_discrete(labels = n_lab(dat)) +
  labs(
    title = paste(
      "REM16 and FT1 expression in female (XX) and male (XY) Cannabis sativa"
    ),
    subtitle = "Each point is one plant; panels use independent y-axes",
    x = "Plant sex (chromosomal)",
    y = Y_LABEL
  ) +
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(FIGURES_DIR, "Fig3_REM16_FT1_combined.png"),
  plot = p3,
  width = 8.5,
  height = 5.5,
  dpi = 300
)

# ----------------------- 9. Session information ----------------------

writeLines(
  capture.output(sessionInfo()),
  file.path(DATA_DIR, "sessionInfo.txt")
)

# ----------------------- 10. Finish ----------------------------------

cat("\nAnalysis completed successfully.\n")
cat("CSV/TXT outputs:", normalizePath(DATA_DIR), "\n")
cat("PNG figures:", normalizePath(FIGURES_DIR), "\n")
