
# new_chain="B"
# parameter="1seqs_50diff_05temp" 

# for cycle in "1AA" "1BB" "1CC" "1DD"; do
# Rscript "/Users/keonapang/Desktop/conversion.R" $cycle $new_chain $parameter
# done

# for cycle in "1A" "1B" "1C" "1D"; do
# Rscript "/Users/keonapang/Desktop/conversion.R" $cycle $new_chain $parameter
# done

suppressMessages(library(data.table))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))

args <- commandArgs(trailingOnly = TRUE)
cycle <- args[1]
new_chain <- args[2]
parameter <- args[3]

cat("====================================================\n")

# Predicted binder structure 
DIR_IN <- paste0("/Users/keonapang/Desktop/3_ProteinPMNN_original") # <------- NEW 
cat("DIR_IN: ", DIR_IN, "\n")
setwd(DIR_IN)

# OUTPUT FILE 
DIR_OUT <- "/Users/keonapang/Desktop/5_pdb_original"
if (!dir.exists(DIR_OUT)) {dir.create(DIR_OUT, recursive = TRUE)}
outfile <- paste0(DIR_OUT, "/cycle", cycle, "_", parameter,".pdb")


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

orgfile <- paste0("/Users/keonapang/Desktop/1_AlphaFold/pep",cycle,".pdb") 
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

