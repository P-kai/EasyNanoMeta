import pandas as pd
import argparse
import os

def count_matches(df, conditions):
    counts = {}
    for key, args in conditions.items():
        if len(args) == 2:
            search_in, search_for = args
            counts[key] = df[search_in].str.contains(search_for, na=False).sum()
        elif len(args) == 4:
            column_b, search_b, column_e, search_e = args
            counts[key] = (df[column_b].str.contains(search_b, na=False) & df[column_e].str.contains(search_e, na=False)).sum()
        else:
            raise ValueError("Invalid number of arguments in conditions.")
    return counts

def parse_arguments():
    parser = argparse.ArgumentParser(description='Based on ZymoBIOMICS data, calculating the accuracy of kraken2 or centrifuge classification on species and genus level.')
    parser.add_argument('-i', '--file_path', type=str, required=True, help='The path to the CSV file to be processed.')
    parser.add_argument('-o', '--output_path', type=str, required=True, help='Output directory path for saving the results.')
    parser.add_argument('-p', '--prefix', type=str, default='result', help='Prefix for the output file name.')
    return parser.parse_args()

def main():
    args = parse_arguments()

    # 读取文件并手动添加列名
    df = pd.read_csv(args.file_path, sep='\t', header=None)
    
    # 添加列名
    df.columns = ['A', 'B', 'C', 'D', 'E']

    # 定义检查条件
    conditions = {
        "A_LF": ('B', 'Lactobacillus_fermentum'),
        "B_LF": ('E', 'Limosilactobacillus'),
        "C_LF": ('E', 'Limosilactobacillus fermentum'),
        "D_LF": ('B', 'Lactobacillus_fermentum', 'E', 'Limosilactobacillus'),
        "E_LF": ('B', 'Lactobacillus_fermentum', 'E', 'Limosilactobacillus fermentum'),
        "A_BS": ('B', 'Bacillus_subtilis'),
        "B_BS": ('E', 'Bacillus'),
        "C_BS": ('E', 'Bacillus subtilis'),
        "D_BS": ('B', 'Bacillus_subtilis', 'E', 'Bacillus'),
        "E_BS": ('B', 'Bacillus_subtilis', 'E', 'Bacillus subtilis'),
        "A_SA": ('B', 'Staphylococcus_aureus'),
        "B_SA": ('E', 'Staphylococcus'),
        "C_SA": ('E', 'Staphylococcus aureus'),
        "D_SA": ('B', 'Staphylococcus_aureus', 'E', 'Staphylococcus'),
        "E_SA": ('B', 'Staphylococcus_aureus', 'E', 'Staphylococcus aureus'),
        "A_PA": ('B', 'Pseudomonas_aeruginosa'),
        "B_PA": ('E', 'Pseudomonas'),
        "C_PA": ('E', 'Pseudomonas aeruginosa'),
        "D_PA": ('B', 'Pseudomonas_aeruginosa', 'E', 'Pseudomonas'),
        "E_PA": ('B', 'Pseudomonas_aeruginosa', 'E', 'Pseudomonas aeruginosa'),
        "A_SE": ('B', 'Salmonella_enterica'),
        "B_SE": ('E', 'Salmonella'),
        "C_SE": ('E', 'Salmonella enterica'),
        "D_SE": ('B', 'Salmonella_enterica', 'E', 'Salmonella'),
        "E_SE": ('B', 'Salmonella_enterica', 'E', 'Salmonella enterica'),
        "A_EC": ('B', 'Escherichia_coli'),
        "B_EC": ('E', 'Escherichia'),
        "C_EC": ('E', 'Escherichia coli'),
        "D_EC": ('B', 'Escherichia_coli', 'E', 'Escherichia'),
        "E_EC": ('B', 'Escherichia_coli', 'E', 'Escherichia coli'),
        "A_EF": ('B', 'Enterococcus_faecalis'),
        "B_EF": ('E', 'Enterococcus'),
        "C_EF": ('E', 'Enterococcus faecalis'),
        "D_EF": ('B', 'Enterococcus_faecalis', 'E', 'Enterococcus'),
        "E_EF": ('B', 'Enterococcus_faecalis', 'E', 'Enterococcus faecalis'),
        "A_LM": ('B', 'Listeria_monocytogenes'),
        "B_LM": ('E', 'Listeria'),
        "C_LM": ('E', 'Listeria monocytogenes'),
        "D_LM": ('B', 'Listeria_monocytogenes', 'E', 'Listeria'),
        "E_LM": ('B', 'Listeria_monocytogenes', 'E', 'Listeria monocytogenes'),
        "A_CN": ('B', 'Cryptococcus_neoformans'),
        "B_CN": ('E', 'Cryptococcus'),
        "C_CN": ('E', 'Cryptococcus neoformans'),
        "D_CN": ('B', 'Cryptococcus_neoformans', 'E', 'Cryptococcus'),
        "E_CN": ('B', 'Cryptococcus_neoformans', 'E', 'Cryptococcus neoformans'),
        "A_SC": ('B', 'Saccharomyces_cerevisiae'),
        "B_SC": ('E', 'Saccharomyces'),
        "C_SC": ('E', 'Saccharomyces cerevisiae'),
        "D_SC": ('B', 'Saccharomyces_cerevisiae', 'E', 'Saccharomyces'),
        "E_SC": ('B', 'Saccharomyces_cerevisiae', 'E', 'Saccharomyces cerevisiae')
    }

    # 计算满足每个条件的行数
    counts = count_matches(df, conditions)

    # 构建输出文件名
    output_file = os.path.join(args.output_path, f"{args.prefix}_output.csv")

    # 保存结果到 CSV 文件
    df_result = pd.DataFrame.from_dict(counts, orient='index', columns=['Count']).reset_index()

    # 将索引列命名为 'Category'
    df_result.rename(columns={'index': 'Category'}, inplace=True)

    # 将结果保存为 CSV 文件
    df_result.to_csv(output_file, index=False, sep='\t')
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    main()
