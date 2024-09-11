import pandas as pd
import argparse

# 函数：读取文件并合并数据
def merge_files(file1, file2, output_file):
    # 读取文件1
    df1 = pd.read_csv(file1, sep="\t")
    # 提取文件1中 ko 的编号
    df1["ko_number"] = df1["KEGG_Annotation"].apply(lambda x: x.split(":")[1])

    # 读取文件2
    df2 = pd.read_csv(file2, sep="\t", header=0, names=["KO", "PathwayL1", "PathwayL2", "Pathway", "KoDescription"])

    # 按照 ko_number 和 KO 字段匹配并合并
    merged_df = pd.merge(df1, df2, left_on="ko_number", right_on="KO", how="left")

    # 选择需要的列并保存结果
    merged_df[["KEGG_Annotation", "PathwayL1", "PathwayL2", "Pathway", "KoDescription", "Abundance"]].to_csv(output_file, sep="\t", index=False)

# 主函数：处理命令行参数
def main():
    parser = argparse.ArgumentParser(description="合并两个文件，并将文件2的指定列加入文件1。")
    parser.add_argument('-f1', '--file1', required=True, help='输入文件1的路径')
    parser.add_argument('-f2', '--file2', required=True, help='输入文件2的路径')
    parser.add_argument('-o', '--output', required=True, help='输出文件的路径')

    args = parser.parse_args()

    # 合并文件
    merge_files(args.file1, args.file2, args.output)

if __name__ == "__main__":
    main()
