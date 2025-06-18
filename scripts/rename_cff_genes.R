#########

# This script tries to update gene names in a CFF file to their main symbols
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
cff_file <- NULL
HGNC_db <- NULL
NCBI_db <- NULL
out_file <- NULL

# Parse arguments in --param=value format
for (arg in args) {
  if (arg == "--debug") {
    debug_mode <- TRUE
  } else if (grepl("^--cff_file=", arg)) {
    cff_file <- sub("^--cff_file=", "", arg)
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
  cat("Usage: Rscript your_script.R --cff_file=FILE --HGNC_db=FILE --NCBI_db=FILE --out_file=FILE [--debug]\n")
  cat("  --cff_file   : Path to CFF file\n")
  cat("  --HGNC_db    : Path to HGNC database file\n")
  cat("  --NCBI_db    : Path to NCBI database file\n")
  cat("  --out_file   : Path to output file\n")
  cat("  --debug (optional) : Enable debug mode\n")
  stop("Error: Missing required arguments. Please provide all required parameters.\n")
}

# Check if all required arguments are provided
if (is.null(cff_file) || is.null(HGNC_db) || is.null(NCBI_db) || is.null(out_file)) {
  print_usage()  # Print usage message and stop the script
}

# Print the arguments if debug mode is enabled
if (debug_mode) {
  cat("Debug mode is ON\n")
  cat("cff_file:", cff_file, "\n")
  cat("HGNC_db:", HGNC_db, "\n")
  cat("NCBI_db:", NCBI_db, "\n")
  cat("out_file:", out_file, "\n")
}

######## Load CFF file ######## 

column_names <- c(
  "chr1",       # Head gene chromosome
  "pos1",       # Head gene breakpoint
  "strand1",    # Head gene strand (optional)
  "chr2",       # Tail gene chromosome
  "pos2",       # Tail gene breakpoint
  "strand2",    # Tail gene strand (optional)
  "library",    # Library type (NA/DNA/RNA, optional)
  "sample_name",# Sample name
  "sample_type",# Sample type (NA/Tumor/Normal, optional)
  "disease",    # Disease name
  "tool",       # Fusion caller name
  "split_cnt",  # Junction-crossing reads (optional, default "-1")
  "span_cnt",   # Spanning reads (optional, default "-1")
  "t_gene1",    # Head gene symbol/alias
  "t_gene_id1",    
  "t_area1",    # Head gene region (optional)
  "t_gene2",    # Tail gene symbol/alias
  "t_gene_id2",    
  "t_area2"     # Tail gene region (optional)
)

cff <- read.delim(cff_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names=column_names)

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
    #if (debug_mode) cat(symbol, "is HGNC main symbol", "\n")
    return(symbol)
  }
  
  # Check in HGNC_dict for alias and chromosome match
  matched_HGNC <- HGNC_dict %>%
    filter(alias == symbol) %>%
    filter(chromosome == chr)
  
  # If found, return the first matching main_symbol
  if (nrow(matched_HGNC) > 0) {
    if (debug_mode) cat("alias found in HGNC:", symbol, "->", matched_HGNC$main_symbol[1], "\n")
    return(matched_HGNC$main_symbol[1])
  }
  
  # Check in NCBI_dict for alias and chromosome match
  matched_NCBI <- NCBI_dict %>%
    filter(alias == symbol) %>%
    filter(chromosome == chr)
  
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

convert_genes <- function(gene_field, chr, tool) {
  
  # Remove arriba brackets, like: CR381670.1(142260),MIR3648(300)
  if (tool =="arriba") {
    gene_field <- gsub("\\(\\d+\\)", "", gene_field)
  }
  # Remove integrate slash-delimited paralogs, like: DKFZp434K1323/LOC100288778
  if (tool =="integrate") {
    gene_field <- gsub("/", ",", gene_field)
  }
  
  symbols <- unlist(strsplit(gene_field, ","))
  processed_symbols <- mapply(symbol_to_main, symbols, chr)
  
  return(paste(unique(processed_symbols), collapse = ","))
}

# convert the gene name columns
cff_conv <- cff %>%
  mutate(
    t_gene1 = mapply(convert_genes, t_gene1, chr1, tool),
    t_gene2 = mapply(convert_genes, t_gene2, chr2, tool),
  )

# Write the output to a CSV file
write.table(cff_conv, file = out_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)