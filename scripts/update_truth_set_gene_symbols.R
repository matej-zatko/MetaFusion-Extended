#########

# This script tries to update gene symbols in a truth set file to their main symbols
# using multiple databases.

#########

if (!requireNamespace("AnnotationDbi", quietly = TRUE))
  BiocManager::install("AnnotationDbi")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr", repos = "https://cloud.r-project.org/")
if (!requireNamespace("stringr", quietly = TRUE))
  install.packages("stringr", repos = "https://cloud.r-project.org/")

suppressMessages(library(AnnotationDbi))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Initialize default values
debug_mode <- FALSE
truth_file <- NULL
HGNC_db <- NULL
NCBI_db <- NULL
out_file <- NULL

# Parse arguments in --param=value format
for (arg in args) {
  if (arg == "--debug") {
    debug_mode <- TRUE
  } else if (grepl("^--truth_file=", arg)) {
    truth_file <- sub("^--truth_file=", "", arg)
  } else if (grepl("^--HGNC_db=", arg)) {
    HGNC_db <- sub("^--HGNC_db=", "", arg)
  } else if (grepl("^--NCBI_db=", arg)) {
    NCBI_db <- sub("^--NCBI_db=", "", arg)
  } else if (grepl("^--out_file=", arg)) {
    out_file <- sub("^--out_file=", "", arg)
  }
}

# Function to print the usage message
print_usage <- function() {
  cat("Usage: Rscript your_script.R --truth_file=FILE --HGNC_db=FILE --NCBI_db=FILE --out_file=FILE [--debug]\n")
  cat("  --truth_file   : Path to CFF file\n")
  cat("  --HGNC_db    : Path to HGNC database file\n")
  cat("  --NCBI_db    : Path to NCBI database file\n")
  cat("  --out_file   : Path to output file\n")
  cat("  --debug (optional) : Enable debug mode\n")
  stop("Error: Missing required arguments. Please provide all required parameters.\n")
}

# Check if all required arguments are provided
if (is.null(truth_file) || is.null(HGNC_db) || is.null(NCBI_db) || is.null(out_file)) {
  print_usage()  # Print usage message and stop the script
}

# Print the arguments if debug mode is enabled
if (debug_mode) {
  cat("Debug mode is ON\n")
  cat("truth_file:", truth_file, "\n")
  cat("HGNC_db:", HGNC_db, "\n")
  cat("NCBI_db:", NCBI_db, "\n")
  cat("out_file:", out_file, "\n")
}


######## Load HGNC database ######## 

# Read the tab-separated file
df <- read.delim(HGNC_db, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

HGNC_main_set <- unique(c(df$symbol))
HGNC_main_set <- HGNC_main_set[HGNC_main_set != ""]

# Expand rows without using tidyr
df <- df %>%
  mutate(chromosome = str_extract(location, "^[0-9XY]+")) %>%
  mutate(all_aliases = paste(prev_symbol, alias_symbol, sep = "|")) %>%
  mutate(all_aliases = strsplit(all_aliases, "\\|")) 

df <- data.frame(
  alias = unlist(df$all_aliases),
  chromosome = rep(df$chromosome, sapply(df$all_aliases, length)),
  main_symbol = rep(df$symbol, sapply(df$all_aliases, length)),
  stringsAsFactors = FALSE
) %>%
  filter(alias != "") %>%
  unique() %>%
  arrange(main_symbol)

# Final result
HGNC_dict <- df

######## Load NCBI database ######## 

# Read the tab-separated file
df <- read.delim(NCBI_db, sep = "\t", header = TRUE, stringsAsFactors = FALSE)  %>%
  mutate(all_aliases = strsplit(Synonyms, "\\|")) 

NCBI_dict <- data.frame(
  alias = unlist(df$all_aliases),
  chromosome = rep(df$chromosome, sapply(df$all_aliases, length)),
  main_symbol = rep(df$Symbol, sapply(df$all_aliases, length)),
  stringsAsFactors = FALSE
) %>%
  filter(alias != "-") %>%
  unique() %>%
  arrange(main_symbol)

######## Load org.Hs.eg.db database ######## 

df <- AnnotationDbi::select(org.Hs.eg.db, 
                            keys = keys(org.Hs.eg.db, keytype = "ALIAS"), 
                            columns = c("ALIAS", "SYMBOL"), 
                            keytype = "ALIAS")

OrgDB_dict <- df %>%
  rename(alias = ALIAS) %>%
  rename(main_symbol = SYMBOL) %>%
  unique() %>%
  arrange(main_symbol)

symbol_to_main <- function(symbol, chr) {
  # First check if symbol exists in HGNC_dict
  if (symbol %in% HGNC_main_set) {
    if (debug_mode) cat(symbol, "is HGNC main symbol", "\n")
    return(symbol)
  }
  
  # Check in HGNC_dict for alias and chromosome match
  matched_HGNC <- HGNC_dict %>%
    filter(alias == symbol) %>%
    filter(if (!is.na(chr)) chromosome == chr else TRUE)
  
  # If found, return the first matching main_symbol
  if (nrow(matched_HGNC) > 0) {
    if (debug_mode) cat("alias found in HGNC:", symbol, "->", matched_HGNC$main_symbol[1], "\n")
    return(matched_HGNC$main_symbol[1])
  }
  
  # Check in NCBI_dict for alias and chromosome match
  matched_NCBI <- NCBI_dict %>%
    filter(alias == symbol) %>%
    filter(if (!is.na(chr)) chromosome == chr else TRUE)
  
  # If found, return the first matching main_symbol
  if (nrow(matched_NCBI) > 0) {
    if (debug_mode) cat("alias found in NCBI:", symbol, "->", matched_NCBI$main_symbol[1], "\n")
    return(matched_NCBI$main_symbol[1])
  }
  
  matched_OrgDB <- OrgDB_dict %>%
    filter(alias == symbol)
  
  if (nrow(matched_OrgDB) > 0) {
    if (debug_mode) cat("alias found in OrgDB:", symbol, "->", matched_OrgDB$main_symbol[1], "\n")
    return(matched_OrgDB$main_symbol[1])
  }
  
  # If not found, return the input symbol
  if (debug_mode) cat(symbol, "not found", "\n")
  return(symbol)
}

# Read the file
data <- readLines(truth_file)

# Process each line
processed_data <- sapply(data, function(line) {
  parts <- strsplit(line, "\\|")[[1]]
  sample <- parts[1]
  genes <- strsplit(parts[2], "--")[[1]]
  converted_genes <- mapply(symbol_to_main, genes, NA)
  return(paste0(sample, "|", paste(converted_genes, collapse = "--")))
})

# Save the output
writeLines(processed_data, out_file)

cat("Processing complete.\n")