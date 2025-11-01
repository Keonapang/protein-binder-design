# Convert json hotspots to .txt inpu file

# Load required package
library(jsonlite)

# Define the input and output directories
input_dir <- "/mnt/nfs/rigenesmb/biobanks/UKBIOBANK/pangk/Keona_scripts/generative-protein-binder-design/1TNF_hotspots_jsons"
output_file <- "/mnt/nfs/rigenesmb/biobanks/UKBIOBANK/pangk/Keona_scripts/generative-protein-binder-design/input/target_file_1TNF_surface.txt"

# List all JSON files in the directory
json_files <- list.files(input_dir, pattern = "\\.json$", full.names = TRUE)

# Initialize an empty vector to store output lines
output_lines <- c()

# Loop through each JSON
for (f in json_files) {
  # Read and parse JSON
  json_data <- fromJSON(f)
  
  # Extract chain names (like "A"); usually just one
  chains <- names(json_data$hotspot_residues)
  
  # For each chain, build residue strings (e.g., A6, A7, ...)
  for (chain in chains) {
    residues <- json_data$hotspot_residues[[chain]]
    
    # Combine into format [A6,A7,...]
    residue_labels <- paste0(chain, residues)
    line <- paste0("[", paste(residue_labels, collapse = ","), "]")
    
    # Append to output
    output_lines <- c(output_lines, line)
  }
}

# Write all lines to the output file
writeLines(output_lines, output_file)

# Print confirmation
cat("✅ Combined", length(output_lines), "entries into:\n", output_file, "\n")#