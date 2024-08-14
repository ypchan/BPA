#!/usr/bin/env Rscript
# Load commandArgs library for command line argument parsing

# Define the function to manipulate the tree string
mid_root_iqtree <- function(tree_string, outgroup_label) {
    pattern <- paste("(", outgroup_label, ":\\d+\\.\\d+)", sep = "")
    matches <- regmatches(tree_string, gregexpr(pattern, tree_string))
    
    # Extracted string
    extracted_string <- matches[[1]]
    
    # Remove the extracted portion from the input string
    modified_string <- gsub(paste0(',', extracted_string), "", tree_string)
    
    # Remove trailing semicolon if it exists
    modified_string <- sub(";$", "", modified_string)
    
    # Extract the numeric value (c) and divide it by 2
    outgroup_branch_len <- as.numeric(strsplit(extracted_string, ':')[[1]][2])
    outgroup_branch_len_0.5 <- outgroup_branch_len / 2
    
    # Create the final formatted string
    final_string <- paste("(", modified_string, ":", outgroup_branch_len_0.5, ",", outgroup_label, ":", outgroup_branch_len_0.5, ");", sep = "")
    
    return(final_string)
}

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check the number of arguments
if (length(args) != 6) {
  cat("Usage: Rscript script.R --tree iqtree_ml.treefile --outgroup outgroup_label --out iqtree_ml.binary.treefile\n")
  quit(status = 1)
}

# Parse command line options
tree_file_index <- which(args == "--tree")
outgroup_index <- which(args == "--outgroup")
out_file_index <- which(args == "--out")

# Check if all required options are provided
if (length(tree_file_index) == 0 || length(outgroup_index) == 0 || length(out_file_index) == 0) {
  cat("Missing required options.\n")
  quit(status = 1)
}

# Get the file paths and outgroup label
tree_file <- args[tree_file_index + 1]
outgroup_label <- args[outgroup_index + 1]
out_file <- args[out_file_index + 1]

# Read the content of the tree file
input_tree <- readLines(tree_file, warn = FALSE)
if (length(input_tree) == 0) {
  cat("Failed to read the tree file.\n")
  quit(status = 1)
}

# Call the function to manipulate the tree string
result_string <- mid_root_iqtree(input_tree, outgroup_label)

# Write the result to the output file
cat(result_string, file = out_file)
