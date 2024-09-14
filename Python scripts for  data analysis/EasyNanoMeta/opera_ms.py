import os
import subprocess
import argparse

def create_base_output_directory():
    """
    在当前路径下创建名为 'easynanometa_out' 的文件夹
    """
    base_dir = os.path.join(os.getcwd(), "easynanometa_out")
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
        print(f"Created directory: {base_dir}")
    else:
        print(f"Directory already exists: {base_dir}")
    return base_dir

def get_fastq_filenames(long_read_folder):
    """
    获取指定文件夹中的所有 .fastq 或 .fq 文件的文件名及其前缀名。

    参数:
    - long_read_folder: 长读数据文件夹的绝对路径

    返回:
    - fastq_files: 包含所有 .fastq 或 .fq 文件名的列表
    """
    fastq_files = [file_name for file_name in os.listdir(long_read_folder) if file_name.endswith((".fastq", ".fq"))]
    return fastq_files

def opera_ms(long_read_folder, short_read_folder, num_processors):
    """
    运行 OPERA-MS 软件的命令，使用 Singularity 容器来执行，并处理长读数据文件夹中的每个文件。

    参数:
    - long_read_folder: 长读数据文件夹的绝对路径
    - short_read_folder: 短读数据文件夹的绝对路径
    - num_processors: 使用的处理器数量
    """
    # 获取长读数据文件夹中的所有 .fastq 或 .fq 文件
    fastq_files = get_fastq_filenames(long_read_folder)
    
    if fastq_files:
        for fastq_file in fastq_files:
            # 提取文件的前缀名（去掉扩展名）
            prefix = os.path.splitext(fastq_file)[0]
            long_read_file = os.path.join(long_read_folder, fastq_file)
            
            # 创建输出目录
            output_dir = os.path.join(os.getcwd(), 'easynanometa_out', 'operams_out', prefix)
            os.makedirs(output_dir, exist_ok=True)
            
            # 构建命令
            command = (
                f"singularity exec easynanometa1_2.sif /bin/bash -c "
                f"\"source /tools/Miniconda3/bin/activate operams && "
                f"perl /tools/Software/OPERA-MS/OPERA-MS.pl "
                f"--short-read1 {short_read_folder}/{prefix}*1* "
                f"--short-read2 {short_read_folder}/{prefix}*2* "
                f"--long-read {long_read_file} "
                f"--no-ref-clustering "
                f"--num-processors {num_processors} "
                f"--out-dir {output_dir}\""
            )

            # 执行命令
            subprocess.run(command, shell=True, check=True)
            print(f"OPERA-MS command executed with long-read file: {long_read_file} and output directory: {output_dir}")
    else:
        print("No fastq/fq files found in the long read folder.")

def opera(args):
    """
    执行 create_base_output_directory, get_fastq_filenames 和 opera_ms 函数
    """
    # 创建基本输出目录
    create_base_output_directory()

    # 运行 opera_ms 函数
    opera_ms(args.long_read_folder, args.short_read_folder, args.num_processors)

def main():
    """
    主函数，解析命令行参数并调用 opera 函数。
    """
    parser = argparse.ArgumentParser(description="Run OPERA-MS with given long and short read folders.")
    subparsers = parser.add_subparsers(dest='command')

    # opera 子命令
    opera_parser = subparsers.add_parser('opera', help='Run OPERA-MS')
    opera_parser.add_argument('-lf', '--long-read-folder', required=True, help="Absolute path to the long read data folder.")
    opera_parser.add_argument('-sf', '--short-read-folder', required=True, help="Absolute path to the short read data folder.")
    opera_parser.add_argument('-t', '--num-processors', type=int, default=24, help="Number of processors to use (default: 24).")
    
    args = parser.parse_args()
    
    if args.command == 'opera':
        opera(args)

if __name__ == "__main__":
    main()
