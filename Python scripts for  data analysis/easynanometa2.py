import os
import subprocess
import argparse

# Flye功能
def run_flye(input_dir, output_dir, singularity_image, flye_executable, threads=24):
    # 检查输入目录是否存在
    if not os.path.exists(input_dir):
        raise ValueError(f"输入目录 {input_dir} 不存在")
    
    # 检查输出目录是否存在，如果不存在则创建
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 遍历输入目录中的所有fastq文件
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.fastq'):
            input_file = os.path.join(input_dir, file_name)
            sample_name = os.path.splitext(file_name)[0]  # 获取文件名（去掉扩展名）
            output_path = os.path.join(output_dir, f"{sample_name}_flye")

            # 如果输出目录不存在，则创建
            if not os.path.exists(output_path):
                os.makedirs(output_path)

            # 生成 Singularity Flye 命令
            command = [
                "singularity", "exec", "-B",
                f"{input_dir}:/mnt",
                singularity_image,
                flye_executable,
                "--nano-raw", f"/mnt/{file_name}",
                "--threads", str(threads),
                "--meta",
                "-o", f"/mnt/{sample_name}_flye"
            ]

            # 打印命令以检查是否正确
            print(f"Running Flye on {file_name} with output to {output_path}")
            print("Command:", " ".join(command))

            # 执行命令
            try:
                subprocess.run(command, check=True)
                print(f"Flye assembly for {file_name} completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error running Flye on {file_name}: {e}")

# Kraken2功能
def run_kraken2(input_dir, output_dir, singularity_image, kraken2_executable, kraken2_db, threads=24):
    # 检查输入目录是否存在
    if not os.path.exists(input_dir):
        raise ValueError(f"输入目录 {input_dir} 不存在")
    
    # 检查输出目录是否存在，如果不存在则创建
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 遍历输入目录中的所有fastq文件
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.fastq'):
            input_file = os.path.join(input_dir, file_name)
            sample_name = os.path.splitext(file_name)[0]  # 获取文件名（去掉扩展名）

            # 生成输出文件路径
            report_file = os.path.join(output_dir, f"{sample_name}_kraken2_report")
            result_file = os.path.join(output_dir, f"{sample_name}_kraken2_result")

            # 生成 Singularity Kraken2 命令
            command = [
                "singularity", "exec", "-B",
                f"{kraken2_db}:/mnt", singularity_image,
                "/bin/bash", "-c",
                f"source /tools/Miniconda3/bin/activate base && "
                f"{kraken2_executable} --db /mnt/k2_standard/ "
                f"--report {report_file} "
                f"--output {result_file} "
                f"-t {threads} {input_file}"
            ]

            # 打印命令以检查是否正确
            print(f"Running Kraken2 on {file_name} with report output to {report_file}")
            print("Command:", " ".join(command))

            # 执行命令
            try:
                subprocess.run(command, check=True)
                print(f"Kraken2 analysis for {file_name} completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error running Kraken2 on {file_name}: {e}")

# 主函数，选择Flye或Kraken2
def main():
    # 创建顶层解析器，仅显示flye和kraken2的选项
    parser = argparse.ArgumentParser(prog='easynano', description='Batch execute Flye or Kraken2 for multiple fastq files.')
    subparsers = parser.add_subparsers(dest='tool', help='Choose the tool to run.')

    # 创建flye的子命令解析器
    parser_flye = subparsers.add_parser('flye', help='Run Flye assembler')
    parser_flye.add_argument('-i', '--input_dir', required=True, help='Directory containing the input fastq files.')
    parser_flye.add_argument('-o', '--output_dir', required=True, help='Directory for storing the Flye output.')
    parser_flye.add_argument('-s', '--singularity_image', required=True, help='Path to the Singularity image (easynanometa1_2.sif).')
    parser_flye.add_argument('-f', '--flye_executable', default='/tools/Software/Flye-2.9.2/bin/flye', help='Path to Flye executable within the Singularity image.')
    parser_flye.add_argument('-t', '--threads', default=24, type=int, help='Number of threads to use.')

    # 创建kraken2的子命令解析器
    parser_kraken2 = subparsers.add_parser('kraken2', help='Run Kraken2 classifier')
    parser_kraken2.add_argument('-i', '--input_dir', required=True, help='Directory containing the input fastq files.')
    parser_kraken2.add_argument('-o', '--output_dir', required=True, help='Directory for storing the Kraken2 output.')
    parser_kraken2.add_argument('-s', '--singularity_image', required=True, help='Path to the Singularity image (easynanometa1_2.sif).')
    parser_kraken2.add_argument('-k', '--kraken2_executable', default='/tools/Software/kraken2-2.1.3/kraken2', help='Path to Kraken2 executable within the Singularity image.')
    parser_kraken2.add_argument('-db', '--kraken2_db', default='/backup/database/kraken2', help='Path to the Kraken2 database.')
    parser_kraken2.add_argument('-t', '--threads', default=24, type=int, help='Number of threads to use.')

    args = parser.parse_args()

    # 选择工具
    if args.tool == 'flye':
        run_flye(args.input_dir, args.output_dir, args.singularity_image, args.flye_executable, args.threads)
    elif args.tool == 'kraken2':
        run_kraken2(args.input_dir, args.output_dir, args.singularity_image, args.kraken2_executable, args.kraken2_db, args.threads)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
