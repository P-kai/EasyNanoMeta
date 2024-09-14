import os
import argparse
import subprocess
import sys

def kraken2(kraken2_db, num_threads):
    # 创建kraken2文件夹
    kraken2_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "kraken2_out")
    if not os.path.exists(kraken2_out_folder):
        os.makedirs(kraken2_out_folder)
        print(f"Created folder: {kraken2_out_folder}")
    else:
        print(f"Folder already exists: {kraken2_out_folder}")
    
    # 固定 host_removal_out_folder 路径为 ./easynanometa_result/host_removal_out
    host_removal_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "host_removal_out")
    
    # 遍历 host_removal_out 文件夹中的所有_fitted_raw.fastq文件
    for fastq_file in os.listdir(host_removal_out_folder):
        if fastq_file.endswith("_fitted_raw.fastq"):
            file_prefix = os.path.splitext(fastq_file)[0]
            fastq_path = os.path.join(host_removal_out_folder, fastq_file)
            # 生成输出的文件路径
            kraken2_report_path = os.path.join(kraken2_out_folder, f"{file_prefix.replace('_fitted_raw', '')}_kraken2_report")
            kraken2_classification_path = os.path.join(kraken2_out_folder, f"{file_prefix.replace('_fitted_raw', '')}_kraken2_result")
                   
            # 使用kraken2进行三代宏基因组物种注释
            kraken2_cmd = (
                f"singularity exec -B {kraken2_db}:/mnt easynanometa1_8.sif /bin/bash -c \"source /tools/Miniconda3/bin/activate base && /tools/Software/kraken2-2.1.3/kraken2 "
                f"--db /mnt "
                f"--report {kraken2_report_path} "
                f"--output {kraken2_classification_path} "
                f"--threads {num_threads} {fastq_path} \""
            )   
            # 执行kraken2命令
            try:
                print(f"Running Kraken2 for {fastq_file}...")
                subprocess.run(kraken2_cmd, shell=True, check=True)
                print(f"Kraken2 classification completed for {fastq_file}")
            except subprocess.CalledProcessError as e:
                print(f"Error during Kraken2 classification for {fastq_file}: {e}")
                continue

# 主函数
def main():
    parser = argparse.ArgumentParser(
        description="Find .fastq or .fq files in the specified folder and process them."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=24,
        help="Number of threads to use for all operations (default: 24)."
    )
    parser.add_argument(
        "-kraken2-db", "--kraken2-db",
        required=True,
        help="Path to Kraken2 database path."
    )
    
    # 解析参数
    args = parser.parse_args()
    
    
    # 如果没有提供任何参数（即仅输入了脚本名），打印帮助信息
    if not any(arg.startswith('-') for arg in sys.argv[1:]):
        parser.print_help()
        sys.exit(1)
    
    # 进行Kraken2分析
    kraken2(args.kraken2_db, args.threads)

if __name__ == "__main__":
    main()
    
