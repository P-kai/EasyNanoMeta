import pandas as pd
import argparse

def merge_files(input_file1, input_file2, output_file):
    # Read the first file (main file), assuming it's tab-separated
    df1 = pd.read_csv(input_file1, sep='\t', header=None, engine='python')

    # Read the second file, assuming it's space-separated, using a raw string to avoid warning
    df2 = pd.read_csv(input_file2, sep=r'\s+', header=None, engine='python')

    # Extract the prefix from the second column in file1 up to the second "_" and insert as a new column
    df1.insert(0, 'New_Column', df1[0].str.extract(r'(contig_\d+)'))

    # Ensure the second file has two columns, the first column being contig values and the second being numeric values
    if df2.shape[1] < 2:
        raise ValueError("Input file 2 must contain two columns: contig and value.")

    # Print column names to ensure they are correct
    print("df1 columns:", df1.columns)
    print("df2 columns:", df2.columns)

    # Merge the contents of the second file into the first file based on the "New_Column" column
    merged_df = pd.merge(df1, df2, left_on='New_Column', right_on=df2.columns[0], how='left')

    # Print the merged DataFrame's column names
    print("Merged DataFrame columns:", merged_df.columns)

    # Drop the redundant '_y' column and rename the '1_y' column to 'Value_from_file2'
    merged_df = merged_df.drop(columns=['0_y']).rename(columns={'1_y': 'Value_from_file2'})

    # Keep the 2nd column (from the original file), the 3rd column, and the 14th column (Value_from_file2)
    final_df = merged_df.iloc[:, [1, 2, 13]]

    # Set column names to CDS_ID, Protein_ID, and Abundance
    final_df.columns = ['CDS_ID', 'Protein_ID', 'Abundance']

    # Remove ".1" from Protein_ID column
    final_df['Protein_ID'] = final_df['Protein_ID'].str.replace(r'\.1$', '', regex=True)

    # Save the merged result into a new txt file with tab separation, including column headers
    final_df.to_csv(output_file, sep='\t', header=True, index=False)
    print(f"Merged file saved to {output_file}")

if __name__ == "__main__":
    # Set up command line argument parser
    parser = argparse.ArgumentParser(description='Merge two files based on contig values.')
    parser.add_argument('-i1', '--input_file1', required=True, help='Path to the first input file (main file).')
    parser.add_argument('-i2', '--input_file2', required=True, help='Path to the second input file (contig and values).')
    parser.add_argument('-o', '--output_file', required=True, help='Path to the output txt file.')

    args = parser.parse_args()

    # Merge the files
    merge_files(args.input_file1, args.input_file2, args.output_file)
