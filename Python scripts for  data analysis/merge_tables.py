import pandas as pd
import os
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir', '-i', type=str, default='.', help='Input data directory.')
    parser.add_argument('--identifier', '-d', type=str, default='_abundance', help='The uniform identifier for name of the prepared files.')
    parser.add_argument('--prefix', '-p', type=str, default='result', help='The prefix of output result.')
    parser.add_argument('--output', '-o', type=str, default='.', help='Output data directory.')
    parser.add_argument('--column', '-c', type=str, default='Gene', help='Column name used as the key for merging between Dataframes.')
    opt = parser.parse_args()

    input_dir = opt.input_dir
    identifier = opt.identifier
    mark_column = opt.column
    file_paths = [os.path.join(input_dir, file_name) for file_name in os.listdir(input_dir) if identifier in file_name]

    if not file_paths:
        print(f"No files found with identifier '{identifier}' in directory '{input_dir}'")
        return

    # Read the first file to initialize the merged dataframe
    try:
        merged_df = pd.read_table(file_paths[0], sep="\t", header=0)
        print(f"Initialized with file: {file_paths[0]}")
    except FileNotFoundError as e:
        print(f"Error reading file: {e}")
        return

    # Merge the rest of the files
    for file_path in file_paths[1:]:
        try:
            df = pd.read_table(file_path, sep="\t", header=0)
            print(f"Merging file: {file_path}")
            merged_df = merged_df.merge(df, on=mark_column, how="outer", suffixes=('', f'_{os.path.basename(file_path).split(".")[0]}'))
        except FileNotFoundError as e:
            print(f"Error reading file: {e}")
            continue

    merged_df.fillna(0, inplace=True)
    pre = opt.prefix
    out_dir = opt.output
    output_path = os.path.join(out_dir, f'{pre}_combine.tsv')
    merged_df.to_csv(output_path, sep='\t', index=False)
    print(f"Output saved to {output_path}")

if __name__ == '__main__':
    main()
