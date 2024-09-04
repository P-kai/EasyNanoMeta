import re
import argparse

def process_file(input_file, output_file):
    processed_data = []

    # 读取输入文件
    with open(input_file, 'r') as file:
        for line in file:
            match = re.search(r'^>(\S+).*?(\[.*?\])', line)
            if match:
                accession = match.group(1)
                species = match.group(2).replace(' ', '_')  # 替换species中的空格为下划线
                processed_data.append(f"{accession} {species}")

    # 将处理结果写入输出文件
    with open(output_file, 'w') as file:
        for entry in processed_data:
            file.write(entry + "\n")

if __name__ == "__main__":
    # 设置命令行参数
    parser = argparse.ArgumentParser(description='Process a file to extract accession numbers and species with brackets, replacing spaces with underscores in species.')
    parser.add_argument('-i', '--input', required=True, help='Path to the input file.')
    parser.add_argument('-o', '--output', required=True, help='Path to the output file.')

    args = parser.parse_args()

    # 处理文件
    process_file(args.input, args.output)
