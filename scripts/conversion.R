
# cycle="1A"
# diffusion="50"
# temp="0.5"
# DIR_WORK="/Users/keonapang/Desktop/NVIDIA/Sept12"

# for cycle in "1A" "1B" "1C" "1D"; do
# Rscript "/Users/keonapang/Desktop/conversion.R" $cycle $diffusion $temp $DIR_WORK
# done

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

args <- commandArgs(trailingOnly = TRUE)
cycle <- args[1]
diffusion <- args[2]
temp <- args[3]
DIR_WORK <- args[4]

new_chain <- "B"
cat("====================================================\n")

# Predicted binder structure 
DIR_IN <- DIR_WORK
cat("DIR_IN: ", DIR_IN, "\n")
setwd(DIR_IN)

# OUTPUT FILE 
DIR_OUT <- paste0(DIR_WORK,"/prodigy_input_pdb")
if (!dir.exists(DIR_OUT)) {dir.create(DIR_OUT, recursive = TRUE)}
outfile <- paste0(DIR_OUT, "/cycle", cycle, "_", diffusion, "diff_", temp, "temp.pdb")

# Search for the specific file
matching_file <- list.files(
  path = DIR_IN,  
  pattern = paste0(".*", cycle,".*rfdiffusion*\\.pdb$"), # <------- NEW 
  full.names = TRUE
)

if (length(matching_file) > 0) {
  infile <- matching_file[1] 
} else {
  infile <- NULL
  cat("No matching file found.\n")
}

# Target ApoB AlphaFold2 structure (.pdb)

orgfile <- paste0(DIR_WORK, "/pep",cycle,".pdb") 
cat("orgfile: ", orgfile, "\n")
cat("infile: ", infile, "\n")
cat("Parameter: ", parameter, "\n")

#########################################################

modify_pdb_chain <- function(new_chain, infile, orgfile) {
  
  # Read the PDB file into R
  full_path <- infile
  pdb_data <- readLines(full_path)
  
  # Initialize an empty vector
  modified_data <- c()
  
  # Loop through each line in the PDB file
  for (line in pdb_data) {
    if (startsWith(line, "ATOM") || startsWith(line, "HETATM")) {
      substr(line, 22, 22) <- new_chain
    }
    modified_data <- c(modified_data, line)
  }
  # remove the first line "MODEL"
  modified_data <- tail(modified_data, -1)
  
  target_seq <- orgfile
  target_data <- readLines(target_seq)
  
  # Remove the last two lines 
  if (length(target_data) > 2) {
    target_data <- head(target_data, -2) 
  }
  
  # Append target_seq + modified_data row-wise 
  combined_data <- c(target_data, modified_data)
  
  # Output file name with "_new" appended
  writeLines(combined_data, outfile)
  cat("COMPLETED:", outfile, "\n")
}

modify_pdb_chain(new_chain, infile, orgfile)
cat("====================================================\n\n")

