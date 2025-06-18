#########

# This script tries to update Ensembl gene IDs in a CFF file
# to their current version.

#########

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr", repos = "https://cloud.r-project.org/")
}
if (!requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr", repos = "https://cloud.r-project.org/")
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", repos = "https://cloud.r-project.org/")
}

suppressMessages(library(httr))
suppressMessages(library(jsonlite))
suppressMessages(library(dplyr))

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Initialize default values
debug_mode <- FALSE
cff_file <- NULL
out_file <- NULL

# Parse arguments in --param=value format
for (arg in args) {
  if (arg == "--debug") {
    debug_mode <- TRUE
  } else if (grepl("^--cff_file=", arg)) {
    cff_file <- sub("^--cff_file=", "", arg)
  } else if (grepl("^--out_file=", arg)) {
    out_file <- sub("^--out_file=", "", arg)
  }
}

# Function to print the usage message
print_usage <- function() {
  cat("Usage: Rscript your_script.R --cff_file=FILE --out_file=FILE [--debug]\n")
  cat("  --cff_file   : Path to CFF file\n")
  cat("  --out_file   : Path to output file\n")
  cat("  --debug (optional) : Enable debug mode\n")
  stop("Error: Missing required arguments. Please provide all required parameters.\n")
}

# Check if all required arguments are provided
if (is.null(cff_file) || is.null(out_file)) {
  print_usage()  # Print usage message and stop the script
}

# Print the arguments if debug mode is enabled
if (debug_mode) {
  cat("Debug mode is ON\n")
  cat("cff_file:", cff_file, "\n")
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

# Check the validity of ID columns in CFF 
cff <- cff %>%
  mutate(
    t_gene_id1 = ifelse(grepl("^ENSG", t_gene_id1), t_gene_id1, NA),  # Replace with NA if it doesn't start with "ENSG"
    t_gene_id2 = ifelse(grepl("^ENSG", t_gene_id2), t_gene_id2, NA)   # Replace with NA if it doesn't start with "ENSG"
  )

# Define Ensembl REST API endpoint
server <- "https://rest.ensembl.org"
ext <- "/archive/id"
ensg_in_cff <- unique(c(cff$t_gene_id1, cff$t_gene_id2))

# Function to process a batch of IDs
process_batch <- function(input_ids) {
  
  # Try to send the request and handle possible failures
  tryCatch({
    r <- POST(
      paste0(server, ext), 
      content_type("application/json"), 
      accept("application/json"), 
      body = toJSON(list(id = input_ids))  # Pass the IDs as JSON
    )
    
    # Check if the request was successful
    stop_for_status(r)
    
    # Convert response content to parsed JSON
    json_data <- content(r, as = "parsed", type = "application/json")
    
    # Extract "id" and all "possible_replacement" stable_ids (comma-separated)
    df <- data.frame(
      id = sapply(json_data, `[[`, "id"),
      possible_replacement = sapply(json_data, function(x) {
        if (!is.null(x$possible_replacement) && length(x$possible_replacement) > 0) {
          x$possible_replacement[[1]]$stable_id  # Extract first stable_id
        } else {
          NA  # Handle cases where there is no replacement
        }
      }),
      stringsAsFactors = FALSE
    )
    
    df <- df[!is.na(df$possible_replacement), ]
    
    return(df)
    
  }, error = function(e) {
    # If there's an error, print the message and return NULL
    message("Error occurred while processing the batch: ", e$message)
    return(NULL)
  })
}

# Function to handle retries
retry_request <- function(input_ids, retries = 3, delay = 5) {
  attempt <- 1
  result <- NULL
  
  while (attempt <= retries && is.null(result)) {
    message("Attempt ", attempt, " of ", retries, "...")
    result <- process_batch(input_ids)
    if (is.null(result)) {
      message("Retrying in ", delay, " seconds...")
      Sys.sleep(delay)  # Wait before retrying
    }
    attempt <- attempt + 1
  }
  
  if (is.null(result)) {
    message("All retry attempts failed. Aborting process entirely.")
  }
  
  return(result)
}

# Initialize an empty list to store results
all_results <- list()

# Process in batches of 1000
batch_size <- 1000
for (i in seq(1, length(ensg_in_cff), by = batch_size)) {
  # Select a batch of unmapped_ensg
  batch_ids <- ensg_in_cff[i:min(i + batch_size - 1, length(ensg_in_cff))]
  
  # Print the number of IDs being processed in this batch
  message("Processing ", length(batch_ids), " IDs in batch ", ceiling(i / batch_size))
  
  # Process the current batch with retry logic
  batch_result <- retry_request(batch_ids)
  
  # If processing is successful, store the result
  if (is.null(batch_result)) {
    # If retry fails for a batch, stop further processing entirely
    message("Process aborted due to failure in batch ", ceiling(i / batch_size))
    break  # Exit the loop and stop the process
  }
  
  # If successful, add the result to the list
  message("Successfully processed batch ", ceiling(i / batch_size))
  all_results[[length(all_results) + 1]] <- batch_result
}

# Combine all batches into a single dataframe (if any batches succeeded)
if (length(all_results) > 0) {
  ensg_mapping <- do.call(rbind, all_results)
} else {
  message("No batches were successfully processed. No changes made to CFF.")
}

if (debug_mode) {
  for (i in 1:nrow(ensg_mapping)) {
    cat(ensg_mapping$id[i], "->", ensg_mapping$possible_replacement[i], "\n")
  }
}

update_gene_id <- function(ensg_id, mapping_df) {
  replacement_map <- setNames(mapping_df$possible_replacement, mapping_df$id)
  ensg_id <- ifelse(ensg_id %in% names(replacement_map), 
                     replacement_map[ensg_id], 
                     ensg_id)
  return(ensg_id)
}


# Use mutate to update t_gene_id1 and t_gene_id2
cff_conv <- cff %>%
  mutate(
    t_gene_id1 = update_gene_id(t_gene_id1, ensg_mapping),
    t_gene_id2 = update_gene_id(t_gene_id2, ensg_mapping)
  )

# Write the updated CFF data to the output file
write.table(cff_conv, file = out_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
