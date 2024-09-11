import pandas as pd
import argparse

def load_species_mapping(mapping_file):
    """Load the ID to Species mapping from a file into a dictionary."""
    species_mapping = {}
    with open(mapping_file, 'r') as f:
        for line in f:
            if line.strip():  # Skip empty lines
                parts = line.split()
                species_mapping[parts[0]] = parts[1]
    return species_mapping

def annotate_diamond_output(diamond_file, species_map, output_file):
    """Annotate DIAMOND output with species information."""
    with open(diamond_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            fields = line.strip().split('\t')
            seq_id = fields[1]  # Assuming the ID is in the second column
            species = species_map.get(seq_id, "[Unknown]")  # Get species from map or use [Unknown]
            outfile.write(f"{line.strip()}\t{species}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Annotate DIAMOND output with species information.')
    parser.add_argument('-d', '--diamond', required=True, help='Path to the DIAMOND output file.')
    parser.add_argument('-m', '--mapping', required=True, help='Path to the ID to Species mapping file.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output file with species information.')

    args = parser.parse_args()

    # Load the ID to Species mapping into a dictionary
    species_map = load_species_mapping(args.mapping)
    
    # Annotate the DIAMOND output with the species information
    annotate_diamond_output(args.diamond, species_map, args.output)
