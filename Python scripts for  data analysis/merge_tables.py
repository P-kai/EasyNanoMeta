import pandas as pd
import os
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_dir','-i', type=str, default='.',help='Input data directory.')
    parser.add_argument('--identifier','-d', type=str, default='_abundance', help='The uniform identifier for name of the prepared files.')
    parser.add_argument('--prefix','-p', type=str, default='result',help='The prefix of output result.')
    parser.add_argument('--output','-o', type=str, default='.', help='Output data directory.')
    parser.add_argument('--column','-c', type=str, default='Gene', help='Column name used as the key for merging between Dataframes.')
    opt = parser.parse_args()

input_dir = opt.input_dir
identifier = opt.identifier
mark_column = opt.column
file_paths = [file_name for file_name in os.listdir(input_dir) if identifier in file_name]
merged_df = pd.DataFrame()
file_path = file_paths[0]
merged_df = pd.read_table(file_path, sep="\t", header=0)
drop_columns = [column + '_y' for column in merged_df.columns]

drop_columns.remove(mark_column + "_y")

for file_path in file_paths:
    df = pd.read_table(file_path, sep="\t", header=0)
    merged_df = merged_df.merge(df, on=mark_column, how="left")

merged_df.fillna(0, inplace=True)
merged_df.drop(columns=drop_columns,inplace=True)
pre = opt.prefix
out_dir = opt.output
merged_df.to_csv(f'{out_dir}/{pre}_combine', sep='\t',index=False)
