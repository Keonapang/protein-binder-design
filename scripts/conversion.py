import os
import sys

# Parse command-line arguments
if len(sys.argv) != 6:
    print("Usage: python script.py <target_file> <binder_file> <diffusion> <temp> <DIR_OUT>")
    sys.exit(1)

target_file = sys.argv[1]
binder_file = sys.argv[2]
diffusion = sys.argv[3]
temp = sys.argv[4]
DIR_OUT = sys.argv[5]

new_chain = "B"
cycle = os.path.splitext(os.path.basename(target_file))[0]

# Add the prefix "5_" to the extracted filename
file_name = os.path.basename(binder_file)  # Extract the file name from the full path
file_part = "_".join(file_name.split("_")[1:])  # Remove everything before and including the first "_"
outname = f"5_{file_part}"

print("==================== ", cycle, "======================")

# OUTPUT FILE
if not os.path.exists(DIR_OUT):
    os.makedirs(DIR_OUT)
outfile = os.path.join(DIR_OUT, outname)

print("target protein PDB: ", target_file)

#########################################################

def modify_pdb_chain(new_chain, binder_file, target_file):
    # Read the PDB file into Python
    with open(binder_file, "r") as f:
        pdb_data = f.readlines()

    # Initialize an empty list
    modified_data = []

    # Loop through each line in the PDB file
    for line in pdb_data:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # Replace the chain identifier (column 22)
            line = line[:21] + new_chain + line[22:]
        modified_data.append(line)

    # Remove the first line "MODEL" (if it exists)
    if modified_data and modified_data[0].strip() == "MODEL":
        modified_data = modified_data[1:]

    # Read the target PDB file
    with open(target_file, "r") as f:
        target_data = f.readlines()

    # Remove the last two lines from the target PDB file
    if len(target_data) > 2:
        target_data = target_data[:-2]

    # Combine target_data and modified_data
    combined_data = target_data + modified_data

    # Write the combined data to the output file
    with open(outfile, "w") as f:
        f.writelines(combined_data)

    print("Merged result:", outfile)

#########################################################

modify_pdb_chain(new_chain, binder_file, target_file)
print("====================================================\n\n")