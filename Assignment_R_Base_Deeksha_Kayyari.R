###############################################################################
# Title: DNA Sequence Analysis (Assignment 2 - Deeksha Kayyari)
# Description:
#   This R program performs the following tasks:
#   1. Reads a DNA sequence from a FASTA file.
#   2. Prints specific nucleotides (10th and 758th).
#   3. Generates and prints reverse complement of the sequence.
#   4. Calculates nucleotide counts per kilobase.
#   5. Converts counts to a structured data frame and analyzes discrepancies.
#
#   Strong error handling is implemented for missing files, invalid input,
#   and incomplete data to ensure robustness and reliability.
###############################################################################

# ============================ LIBRARIES ======================================
if (!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr", repos = "http://cran.us.r-project.org")
}
library(stringr)

###############################################################################
# ============================ PART 1a ========================================
# Read the complete DNA sequence from the file "chr1_GL383518v1_alt.fa".
###############################################################################

read_sequence <- function(file_path) {
  # --- Error handling for file issues ---
  if (!file.exists(file_path)) stop("File not found.")
  if (file.access(file_path, mode = 4) != 0) stop("Permission denied.")
  if (file.info(file_path)$size == 0) stop("The file exists but is empty.")
  
  # --- Read and clean sequence ---
  seq_lines <- readLines(file_path, warn = FALSE)
  seq_lines <- seq_lines[!grepl("^>", seq_lines)]  # remove FASTA header lines
  seq <- paste(seq_lines, collapse = "")
  seq <- toupper(seq)
  seq <- gsub("[^ACGT]", "", seq)  # remove non-ACGT characters
  
  if (nchar(seq) == 0) stop("Sequence is empty after cleaning invalid bases.")
  return(seq)
}

sequence_file <- "chr1_GL383518v1_alt.fa"

sequence <- tryCatch({
  read_sequence(sequence_file)
}, error = function(e) {
  cat("Error reading sequence:", e$message, "\n")
  return(NULL)
})

if (!is.null(sequence)) {
  cat("\n=== PART 1a: Sequence Loaded Successfully ===\n")
  cat("Total sequence length:", nchar(sequence), "bases\n")
}

###############################################################################
# ============================ PART 1b ========================================
# Print the 10th and 758th bases of the DNA sequence.
###############################################################################

if (!is.null(sequence)) {
  cat("\n=== PART 1b: Print Specific Bases ===\n")
  cat("10th base:", substr(sequence, 10, 10), "\n")
  cat("758th base:", substr(sequence, 758, 758), "\n")
}

###############################################################################
# ============================ PART 2a ========================================
# Create the reverse complement of the DNA sequence and print the 79th base.
###############################################################################

reverse_complement <- function(seq) {
  if (is.null(seq)) stop("Sequence not loaded.")
  complement <- chartr("ACGT", "TGCA", seq)
  rev_comp <- paste(rev(strsplit(complement, NULL)[[1]]), collapse = "")
  return(rev_comp)
}

rev_seq <- tryCatch({
  if (!is.null(sequence)) reverse_complement(sequence) else NULL
}, error = function(e) {
  cat("Error creating reverse complement:", e$message, "\n")
  return(NULL)
})

if (!is.null(rev_seq)) {
  cat("\n=== PART 2a: Reverse Complement Created ===\n")
  cat("Reverse complement length:", nchar(rev_seq), "bases\n")
  cat("79th base:", substr(rev_seq, 79, 79), "\n")
}

###############################################################################
# ============================ PART 2b ========================================
# Print the bases from position 500–800 of the reverse complement.
###############################################################################

if (!is.null(rev_seq)) {
  cat("\n=== PART 2b: Reverse Complement Bases 500–800 ===\n")
  range500_800 <- substr(rev_seq, 500, 800)
  cat("Bases from 500–800:\n", range500_800, "\n")
}

###############################################################################
# ============================ PART 3a ========================================
# Create a list containing counts of A, C, G, and T per kilobase.
###############################################################################

count_per_kb <- function(seq, kb_size = 1000) {
  seq_len <- nchar(seq)
  if (seq_len == 0) stop("Empty sequence provided.")
  
  count_dict <- list()
  for (i in seq(1, seq_len, kb_size)) {
    segment <- substr(seq, i, min(i + kb_size - 1, seq_len))
    counts <- list(
      A = str_count(segment, "A"),
      C = str_count(segment, "C"),
      G = str_count(segment, "G"),
      T = str_count(segment, "T")
    )
    count_dict[[as.character(i)]] <- counts
  }
  return(count_dict)
}

counts_dict <- tryCatch({
  if (!is.null(sequence)) count_per_kb(sequence) else NULL
}, error = function(e) {
  cat("Error counting nucleotides:", e$message, "\n")
  return(NULL)
})

if (!is.null(counts_dict)) {
  cat("\n=== PART 3a: Nucleotide Counts per Kilobase ===\n")
  cat("Total kilobases analyzed:", length(counts_dict), "\n")
  cat("\nExample output (first 3 kilobases):\n")
  for (i in 1:min(3, length(counts_dict))) {
    kb <- names(counts_dict)[i]
    cts <- counts_dict[[i]]
    cat("Start", kb, "-> A:", cts$A, "C:", cts$C, "G:", cts$G, "T:", cts$T, "\n")
  }
}

###############################################################################
# ============================ PART 4a ========================================
# Convert the list of counts into a structured data frame.
###############################################################################

convert_to_dataframe <- function(count_list) {
  if (is.null(count_list) || length(count_list) == 0)
    stop("No data to convert.")
  
  df <- data.frame(
    Start_Position = as.integer(names(count_list)),
    A = sapply(count_list, function(x) x$A),
    C = sapply(count_list, function(x) x$C),
    G = sapply(count_list, function(x) x$G),
    T = sapply(count_list, function(x) x$T)
  )
  return(df)
}

nuc_df <- tryCatch({
  if (!is.null(counts_dict)) convert_to_dataframe(counts_dict) else NULL
}, error = function(e) {
  cat("Error creating data frame:", e$message, "\n")
  return(NULL)
})

if (!is.null(nuc_df)) {
  cat("\n=== PART 4a: Data Frame Created ===\n")
  print(head(nuc_df, 5))
}

###############################################################################
# ============================ PART 4b ========================================
# Verify that all kilobases are included.
###############################################################################

if (!is.null(nuc_df)) {
  cat("\n=== PART 4b: Verification of Kilobases ===\n")
  cat("Total rows in data frame:", nrow(nuc_df), "\n")
}

###############################################################################
# ============================ PART 4c ========================================
# Display base count rows (A, C, G, T) for first few kilobases.
###############################################################################

if (!is.null(nuc_df)) {
  cat("\n=== PART 4c: Base Count Rows ===\n")
  print(head(nuc_df[, c("A", "C", "G", "T")], 5))
}

###############################################################################
# ============================ PART 4d ========================================
# Calculate the total bases per kilobase (A + C + G + T).
###############################################################################

if (!is.null(nuc_df)) {
  nuc_df$Sum <- rowSums(nuc_df[, c("A", "C", "G", "T")])
  cat("\n=== PART 4d: Row Sums Calculated ===\n")
  print(head(nuc_df, 5))
}

###############################################################################
# ============================ PART 4e ========================================
# Validate sums and explain discrepancies.
###############################################################################

if (!is.null(nuc_df)) {
  expected_sum <- 1000
  incorrect <- subset(nuc_df, Sum != expected_sum)
  
  cat("\n=== PART 4e: Validation and Explanation ===\n")
  cat("Expected sum per kilobase:", expected_sum, "\n")
  cat("Number of kilobases analyzed:", nrow(nuc_df), "\n")
  
  if (nrow(incorrect) > 0) {
    cat("\nKilobases with unexpected totals:\n")
    print(incorrect)
    cat("\nExplanation:\n")
    cat("1. The final kilobase may contain fewer than 1000 bases if the total sequence\n")
    cat("   length is not divisible by 1000.\n")
    cat("2. Ambiguous bases (N, R, Y, etc.) were removed during cleaning, lowering counts.\n")
    cat("3. Line breaks or FASTA formatting may slightly reduce sequence length.\n")
  } else {
    cat("\nAll kilobases sum to 1000 bases as expected.\n")
  }
}

cat("\nProgram completed successfully.\n")
###############################################################################
