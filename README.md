# Cannabis sativa: REM16 and FT1 Gene Expression Analysis

## Overview

This repository contains an R script for analyzing the expression of two genes in *Cannabis sativa*:

- **REM16:** `LOC115699937`
- **FT1:** `LOC115696989`

The analysis compares gene expression between female (XX) and male (XY) plants using RNA-seq expression data. It generates summary statistics, performs Wilcoxon rank-sum tests, and produces three figures for visualizing gene expression by plant sex.

The script is designed to run directly from the repository without requiring manual modification of file paths.

## 1. Requirements

Before running the analysis, ensure that you have:

- Ubuntu or another compatible Linux environment.
- R and Rscript.
- Git.
- GNU Make and build tools, required if installing R packages (data.table and ggplot2) from source is needed.

### Check whether the required software is installed

Open an Ubuntu terminal and run:

```bash
R --version
Rscript --version
git --version
make --version
```

If a command is not found, follow the installation instructions below.

### Install R and Git

If R or Git is not installed, run:

```bash
sudo apt update
sudo apt install r-base git
```

Verify the installation:

```bash
R --version
Rscript --version
git --version
```

### Install build tools

The R script automatically installs the required R packages if they are missing. However, Ubuntu may need additional build tools to compile packages from source.

If you encounter an error such as:

```text
sh: 1: make: not found
```

install the build tools:

```bash
sudo apt update
sudo apt install build-essential
```

This package provides `make` and essential compilation tools used during package installation.

Verify that `make` is available:

```bash
make --version
```

### Important: Ubuntu password

When Ubuntu asks for a password after you run a command using `sudo`, enter your **Ubuntu user password**.

- It is the password associated with your Ubuntu account.
- It is not your GitHub password.
- It is not necessarily your Windows password.
- While typing, the terminal will not display characters, dots, or asterisks. This is normal. Type the password and press **Enter**.

These system installations generally need to be performed only once per Ubuntu environment.

## 2. Clone the Repository

Clone the repository from GitHub:

```bash
git clone https://github.com/ereenyhozain712/oLLMpiques.git
```

Navigate to the repository directory:

```bash
cd oLLMpiques
```

Check the repository structure:

```bash
ls
```

The relevant files should be organized as follows:

```text
oLLMpiques/
├── code/
│   └── sex_gene_expression.R
├── data/
│   └── 2_Data_RNASeq_Cannabis_Sex.csv
├── figures/
└── README.md
```

The input CSV file must be present in the `data/` directory.

## 3. Run the Analysis

From the repository root, run:

```bash
Rscript code/sex_gene_expression.R
```

The script uses paths relative to the repository, so there is no need to edit file paths manually.

**First run:** Execution may take longer because the required R packages may need to be downloaded and installed.

If package installation fails, consult the troubleshooting section below.

## 4. Input Data

The analysis uses the following input file:

`data/2_Data_RNASeq_Cannabis_Sex.csv`

The expected format is a gene-expression matrix where:

- The first column contains gene or locus identifiers.
- Each remaining column represents a plant sample.
- Each row represents a gene.
- Sample names contain a sex suffix, such as `_XX` or `_XY`.

The script identifies the target genes using their specified locus IDs. Other gene names in the input file do not need to begin with `LOC`.

Samples with unrecognized sex suffixes are excluded from the sex-comparison analysis.

## 5. Analysis Workflow

The script performs the following steps:

1. **Data loading:** Reads the RNA-seq expression matrix.
2. **Gene extraction:** Extracts REM16 and FT1 using their locus identifiers.
3. **Sex classification:** Classifies samples as female (XX) or male (XY) based on their sample names.
4. **Expression quality checks:** Reports missing values, negative values, and expression ranges.
5. **Summary statistics:** Calculates the sample count, median, mean, and standard deviation for each gene and sex group.
6. **Statistical testing:** Performs a Wilcoxon rank-sum test to compare expression between female and male plants for each gene.
7. **Visualization:** Generates three figures showing gene expression by plant sex.
8. **Session recording:** Saves R and package version information to support reproducibility.

## 6. Expected Outputs and Verification

After running the script successfully, it displays the dataset summary, sample counts by sex, and descriptive statistics for REM16 and FT1. It also saves the following output files:

### CSV and TXT files — `data/`

| File | Description |
|---|---|
| `REM16_FT1_extracted.csv` | Expression values for REM16 and FT1, including sample names and sex labels. |
| `REM16_FT1_summary_by_sex.csv` | Summary statistics by gene and sex, including sample count, median, mean, and standard deviation. |
| `sessionInfo.txt` | R version and package information used for the analysis. |

The original input file, `2_Data_RNASeq_Cannabis_Sex.csv`, remains unchanged.

### Figures — `figures/`

| File | Description |
|---|---|
| `Fig1_REM16_by_sex.png` | REM16 expression by sex, with boxplots, individual samples, and the Wilcoxon test p-value. |
| `Fig2_FT1_by_sex.png` | FT1 expression by sex, with boxplots, individual samples, and the Wilcoxon test p-value. |
| `Fig3_REM16_FT1_combined.png` | Combined visualization of both genes, with a separate panel for each gene. |

All figures are saved as PNG files at 300 dpi.

### Verify the output files

Once the script finishes, check that the output files were created:

```bash
ls -lh data
ls -lh figures
```

The CSV and TXT files should appear in `data/`, and the three PNG figures should appear in `figures/`.

## Reproducibility

This repository is designed to support reproducible execution from a fresh GitHub clone. The script uses paths relative to the project directory, installs missing R packages into a user-specific library, saves results in predefined directories, and records session information.

The generated figures and tables can be used to inspect the expression patterns of REM16 and FT1 across female and male *Cannabis sativa* plants.
