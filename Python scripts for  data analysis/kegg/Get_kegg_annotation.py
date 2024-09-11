import pandas as pd
import argparse

def get_kegg_annotation(input_file1, input_file2, output_file):
    # Read file1 and store it in a DataFrame, retaining the header
    df1 = pd.read_csv(input_file1, sep="\t")
    
    # Read file2 and store it in a DataFrame
    df2 = pd.read_csv(input_file2, sep="\t", header=None, names=["Protein_ID", "KEGG_Annotation"])
    
    # Merge the two DataFrames based on "Protein_ID"
    df_merged = pd.merge(df1, df2, on="Protein_ID", how="inner")
    
    # Output the specified columns: CDS_ID, Protein_ID, KEGG_Annotation, Abundance
    df_merged[["CDS_ID", "Protein_ID", "KEGG_Annotation", "Abundance"]].to_csv(output_file, sep="\t", index=False, header=True)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(description="Merge KEGG annotation files and extract results.")
    parser.add_argument('-i1', '--input_file1', required=True, help="Path to the first input file (e.g., 'kegg_test_match_exct.out').")
    parser.add_argument('-i2', '--input_file2', required=True, help="Path to the second input file (e.g., 'Plants.NCBI2KEGG.txt').")
    parser.add_argument('-o', '--output_file', required=True, help="Path to the output file (e.g., 'kplant_test_out.txt').")
    
    args = parser.parse_args()

    # Call the function to process the data
    get_kegg_annotation(args.input_file1, args.input_file2, args.output_file)
