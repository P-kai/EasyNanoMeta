import pandas as pd
import argparse

# Function: Read files and merge data
def merge_files(file1, file2, output_file):
    # Read file 1
    df1 = pd.read_csv(file1, sep="\t")
    # Extract the KO number from file 1
    df1["ko_number"] = df1["KEGG_Annotation"].apply(lambda x: x.split(":")[1])

    # Read file 2
    df2 = pd.read_csv(file2, sep="\t", header=0, names=["KO", "PathwayL1", "PathwayL2", "Pathway", "KoDescription"])

    # Merge based on the ko_number from file1 and KO field from file2
    merged_df = pd.merge(df1, df2, left_on="ko_number", right_on="KO", how="left")

    # Select the required columns and save the results
    merged_df[["KEGG_Annotation", "PathwayL1", "PathwayL2", "Pathway", "KoDescription", "Abundance"]].to_csv(output_file, sep="\t", index=False)

# Main function: Handle command line arguments
def main():
    parser = argparse.ArgumentParser(description="Merge two files and add specific columns from file 2 into file 1.")
    parser.add_argument('-f1', '--file1', required=True, help='Path to the first input file')
    parser.add_argument('-f2', '--file2', required=True, help='Path to the second input file')
    parser.add_argument('-o', '--output', required=True, help='Path to the output file')

    args = parser.parse_args()

    # Merge the files
    merge_files(args.file1, args.file2, args.output)

if __name__ == "__main__":
    main()
