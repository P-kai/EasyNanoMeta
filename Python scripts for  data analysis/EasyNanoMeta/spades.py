import os
import argparse
import subprocess

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

def get_long_fastq_filenames(long_read_folder):
    """
    获取指定文件夹中的所有 .fastq 或 .fq 文件的文件名及其前缀名。

    参数:
    - long_read_folder: 长读数据文件夹的绝对路径

    返回:
    - fastq_files: 包含所有 .fastq 或 .fq 文件名的列表
    """
    fastq_files = [file_name for file_name in os.listdir(long_read_folder) if file_name.endswith((".fastq", ".fq"))]
    return fastq_files

def spades_run(long_read_folder, short_read_folder):
    """
    运行 SPAdes 的函数。根据指定的长读和短读数据文件夹执行命令。

    参数:
    - long_read_folder: 长读数据文件夹的绝对路径
    - short_read_folder: 短读数据文件夹的绝对路径
    """
    # 获取所有长读数据文件
    long_read_files = get_long_fastq_filenames(long_read_folder)
    
    if not long_read_files:
        print("No .fastq or .fq files found in the long read folder.")
        return
    
    # 基础输出目录
    spades_out_dir = os.path.join(os.getcwd(), "easynanometa_out", "spades_out")
    
    if not os.path.exists(spades_out_dir):
        os.makedirs(spades_out_dir)
        print(f"Created base output directory: {spades_out_dir}")
    else:
        print(f"Base output directory already exists: {spades_out_dir}")
    
    # 依次处理每个长读文件
    for long_read_file in long_read_files:
        # 获取前缀名（假设前缀名为文件名的前部分，去掉扩展名）
        prefix = os.path.splitext(long_read_file)[0]
        
        # 为每个前缀创建单独的输出文件夹
        output_dir = os.path.join(spades_out_dir, prefix)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print(f"Created output directory for prefix '{prefix}': {output_dir}")
        
        # 构建命令
        command = (
            f'singularity exec easynanometa1_6.sif /bin/bash -c "'
            f'source /tools/Miniconda3/bin/activate base && '
            f'/tools/Software/SPAdes-4.0.0-Linux/bin/spades.py --meta '
            f'-1 {short_read_folder}/{prefix}*1* '
            f'-2 {short_read_folder}/{prefix}*2* '
            f'--nanopore {long_read_folder}/{long_read_file} '
            f'--meta -t 24 -o {output_dir}"'
        )

        # 执行命令
        print(f"Running command: {command}")
        subprocess.run(command, shell=True)

def spades(args):
    """
    spades 子命令的主函数，依次执行 create_base_output_directory, get_long_fastq_filenames, spades_run 函数。

    参数:
    - args: 包含 long_read_folder 和 short_read_folder 参数的命令行参数对象
    """
    create_base_output_directory()
    spades_run(args.long_read_folder, args.short_read_folder)

def main():
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(description="Run SPAdes with specified long and short read folders.")
    
    # 创建子命令解析器
    subparsers = parser.add_subparsers(title="Subcommands", description="Available subcommands", help="sub-command help")
    
    # 定义 spades 子命令
    parser_spades = subparsers.add_parser("spades", help="Run SPAdes on specified long and short read folders")
    parser_spades.add_argument('-lf', '--long_read_folder', type=str, required=True, help="长读数据文件夹的绝对路径")
    parser_spades.add_argument('-sf', '--short_read_folder', type=str, required=True, help="短读数据文件夹的绝对路径")
    parser_spades.set_defaults(func=spades)
    
    # 解析命令行参数
    args = parser.parse_args()
    
    # 调用子命令对应的函数
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
