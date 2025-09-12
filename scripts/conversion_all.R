
# cycle="1"
# diffusion="50"
# temp="0.5"
# DIR_WORK="/Users/keonapang/Desktop/NVIDIA/Sept12"
# Rscript "/Users/keonapang/Desktop/conversion_all.R" $cycle $diffusion $temp $DIR_WORK

# for cycle in "1A" "1B" "1C" "1D" "2A" "2B" "2C" "2D" "1E" "1F" "1G" "1H" "2E" "2F" "2G" "2H"; do
#     python3.11 /home/ubuntu/3_protein_binder_design.py --cycle "$cycle" --num_seq 1 --diffusion "$diffusion" --temp "$temp"
# done

# for cycle in "1D"; do
#     python3.11 /home/ubuntu/3_protein_binder_design.py --cycle "$cycle" --num_seq 1 --diffusion "$diffusion" --temp "$temp"
# done

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

args <- commandArgs(trailingOnly = TRUE)
cycle <- as.character(args[1])
diffusion <- as.character(args[2])
temp <- as.character(args[3])
DIR_WORK <- as.character(args[4])
# orgfile <- as.character(args[4])

#########################################################
# Define input file names and corresponding chains
#########################################################
setwd(DIR_WORK)

# find matching files in the directory
pattern <- paste0("2_cycle", cycle, "[A-D]_1seqs_", diffusion, "diff_", temp, "temp_rfdiffusion.pdb")
infiles <- list.files(path = ".", pattern = pattern, full.names = TRUE)
if (length(infiles) > 0) {
  print(paste("Found", length(infiles), "file(s):"))
  print(infiles)
} else {
  print("No matching files found.")
}

# if (cycle=="1") {
#   infiles <- c(
#     paste0("2_cycle1A_1seqs_",diffusion,"diff_",temp,"temp_rfdiffusion.pdb"),
#     paste0("2_cycle1B_1seqs_",diffusion,"diff_",temp,"temp_rfdiffusion.pdb"),
#     paste0("2_cycle1C_1seqs_",diffusion,"diff_",temp,"temp_rfdiffusion.pdb"),
#     paste0("2_cycle1D_1seqs_",diffusion,"diff_",temp,"temp_rfdiffusion.pdb")
#   )
# } else {
#   infiles <- c( 
#     paste0("2_cycle2A_1seqs_",diffusion,"diff_",temp,"temp_rfdiffusion.pdb"),
#     paste0("2_cycle2B_1seqs_",diffusion,"diff_",temp,"temp_rfdiffusion.pdb"),
#     paste0("2_cycle2C_1seqs_",diffusion,"diff_",temp,"temp_rfdiffusion.pdb"),
#     paste0("2_cycle2D_1seqs_",diffusion,"diff_",temp,"temp_rfdiffusion.pdb")
#   )
# }

# Corresponding new_chain values for each infile
if (length(infiles)==4){
  new_chains <- c("B", "C", "D", "E")
} else if (length(infiles)==8) {
  new_chains <- c("B", "C", "D", "E", "F", "G", "H", "I")
} else if (length(infiles)==12) {
  new_chains <- c("B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "M", "N")
} else {
  stop("ERROR")
}

# Path to the original file
orgfile <- paste0(DIR_WORK, "/pep", cycle, ".pdb")
cat("\n")
cat("Target structure:", orgfile, "\n\n")

# output directory
DIR_OUT <- paste0(DIR_WORK,"/pymol_pdb")
if (!dir.exists(DIR_OUT)) {dir.create(DIR_OUT, recursive = TRUE)}

# Output file to store the combined result
output_file <- paste0(DIR_OUT, "/cycle", cycle, "_", diffusion, "diff_", temp, "temp_complex.pdb")

#########################################################

# Function to modify the chain and append the data
modify_pdb_chain <- function(new_chain, infile, combined_data) {
  pdb_data <- readLines(infile)
  
  # Modify the chain identifier for lines starting with "ATOM" or "HETATM"
  modified_data <- c()
  for (line in pdb_data) {
    if (startsWith(line, "ATOM") || startsWith(line, "HETATM")) {
      substr(line, 22, 22) <- new_chain
    }
    modified_data <- c(modified_data, line)
  }

  # remove the first line of modified_data that says "MODEL"
    modified_data <- tail(modified_data, -1)  # Keep all lines except the first one

  # Combine the existing data with modified_data
  combined_data <- c(combined_data, modified_data)
  return(combined_data)
}

# Start with the original file as the base data
combined_data <- readLines(orgfile)

# Remove the last two lines from the original file, if they exist
if (length(combined_data) > 2) {
  combined_data <- head(combined_data, -2)  # Keep all lines except the last two
}

# Loop through each infile and new_chain
for (i in seq_along(infiles)) {
  infile <- file.path(DIR_WORK, infiles[i])  # infile <- file.path("/Users/keonapang/Desktop/4_alphafold", infiles[i])

  cat("Processing", infile, "\n\n") 
  new_chain <- new_chains[i]
  
  # Modify the chain and append the infile to combined_data
  combined_data <- modify_pdb_chain(new_chain, infile, combined_data)
}

writeLines(combined_data, output_file)
cat("\n")
cat("Result:", output_file, "\n\n")




