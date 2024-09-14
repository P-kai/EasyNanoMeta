import os
import argparse
import subprocess
import sys

def centrifuge(centrifuge_db, num_threads):
    # 创建centrifuge_out文件夹
    centrifuge_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "centrifuge_out")
    if not os.path.exists(centrifuge_out_folder):
        os.makedirs(centrifuge_out_folder)
        print(f"Created folder: {centrifuge_out_folder}")
    else:
        print(f"Folder already exists: {centrifuge_out_folder}")

    # 固定 host_removal_out_folder 路径为 ./easynanometa_result/host_removal_out
    host_removal_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "host_removal_out")
        
    # 遍历 host_removal_out 文件夹中的所有_fitted_raw.fastq文件
    for fastq_file in os.listdir(host_removal_out_folder):
        if fastq_file.endswith("_fitted_raw.fastq"):
            file_prefix = os.path.splitext(fastq_file)[0]
            fastq_path = os.path.join(host_removal_out_folder, fastq_file)
            # 生成输出的文件路径
            centrifuge_report_path = os.path.join(centrifuge_out_folder, f"{file_prefix.replace('_fitted_raw', '')}_report")
            centrifuge_classification_path = os.path.join(centrifuge_out_folder, f"{file_prefix.replace('_fitted_raw', '')}_result")
                   
            # 步骤1：使用centrifuge进行三代宏基因组物种注释
            centrifuge_db_directory = os.path.dirname(centrifuge_db)
            db_name = os.path.basename(centrifuge_db)
            centrifuge_cmd = (
                f"singularity exec -B {centrifuge_db_directory}:/mnt easynanometa1_7.sif /tools/Software/centrifuge-1.0.4/bin/centrifuge "
                f"-p {num_threads} -x /mnt/{db_name} "
                f"-q {fastq_path} "
                f"--report-file {centrifuge_report_path} "
                f"-S {centrifuge_classification_path}"
            )   
            # 执行centrifuge命令
            try:
                print(f"Running Centrifuge for {fastq_file}...")
                subprocess.run(centrifuge_cmd, shell=True, check=True)
                print(f"Centrifuge classification completed for {fastq_file}")
            except subprocess.CalledProcessError as e:
                print(f"Error during Centrifuge classification for {fastq_file}: {e}")
                continue
    
            # 步骤2：将centrifuge输出结果转换为kraken2结果
            centrifuge_to_kraken_cmd = (
                f"singularity exec -B {centrifuge_db_directory}:/mnt easynanometa1_7.sif /tools/Software/centrifuge-1.0.4/bin/centrifuge-kreport "
                f"-x /mnt/{db_name} "
                f"{centrifuge_classification_path} > "
                f"{centrifuge_out_folder}/{file_prefix}_kraken_report"
            )
            
            # 执行centrifuge-kreport命令
            try:
                print(f"Converting Centrifuge result to Kraken format for {fastq_file}...")
                subprocess.run(centrifuge_to_kraken_cmd, shell=True, check=True)
                print(f"Centrifuge result converted to Kraken report for {fastq_file}")
            except subprocess.CalledProcessError as e:
                print(f"Error during Centrifuge-to-Kraken conversion for {fastq_file}: {e}")
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
        "-centrifuge-db", "--centrifuge-db",
        required=True,
        help="Path to centrifuge database path."
    )
    
    # 解析参数
    args = parser.parse_args()
    
    
    # 如果没有提供任何参数（即仅输入了脚本名），打印帮助信息
    if not any(arg.startswith('-') for arg in sys.argv[1:]):
        parser.print_help()
        sys.exit(1)
    
    # 进行Centrifuge分析
    centrifuge(args.centrifuge_db, args.threads)

if __name__ == "__main__":
    main()
