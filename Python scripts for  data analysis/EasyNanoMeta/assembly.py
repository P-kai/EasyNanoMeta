import os
import subprocess
import argparse
import sys

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

def get_fastq_filenames(directory):
    """
    获取指定目录下所有fastq或fq文件的文件名（包括完整路径）
    """
    filenames = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".fastq") or file.endswith(".fq"):
                # 获取文件的完整路径
                filenames.append(os.path.join(root, file))
    return filenames

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
    
def flye(fastq_file, threads, base_output_dir):
    """
    执行 Flye 命令，使用指定的fastq文件作为输入，并指定线程数
    在指定的 base_output_dir 下创建以文件名为基础的输出文件夹
    """
    filename = os.path.basename(fastq_file)
    file_prefix = os.path.splitext(filename)[0]
    
    flye_out_dir = os.path.join(base_output_dir, "flye_out")
    if not os.path.exists(flye_out_dir):
        os.makedirs(flye_out_dir)
        print(f"Created directory: {flye_out_dir}")
    
    output_dir = os.path.join(flye_out_dir, f"{file_prefix}_flye_out")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")
    else:
        print(f"Directory already exists: {output_dir}")

    command = [
        "singularity", "exec", "easynanometa1_6.sif",
        "/tools/Software/Flye-2.9.2/bin/flye",
        "--nano-raw", fastq_file,
        "--threads", str(threads),
        "--meta",
        "-o", output_dir
    ]
    
    try:
        subprocess.run(command, check=True)
        print(f"Flye command executed successfully for {fastq_file}.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running Flye for {fastq_file}: {e}")

def canu(fastq_file, base_output_dir, genome_size):
    """
    执行 Canu 命令，使用指定的fastq文件作为输入
    在指定的 base_output_dir 下创建以文件名为基础的输出文件夹
    """
    filename = os.path.basename(fastq_file)
    file_prefix = os.path.splitext(filename)[0]
    
    canu_out_dir = os.path.join(base_output_dir, "canu_out")
    if not os.path.exists(canu_out_dir):
        os.makedirs(canu_out_dir)
        print(f"Created directory: {canu_out_dir}")
    
    output_dir = os.path.join(canu_out_dir, file_prefix)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")
    else:
        print(f"Directory already exists: {output_dir}")

    command = [
        "singularity", "exec", "easynanometa1_6.sif",
        "/tools/Software/canu-2.2/bin/canu",
        "-d", output_dir,
        "-p", file_prefix,
        f"genomeSize={genome_size}",
        "-nanopore-raw", fastq_file,
        "minInputCoverage=5"
    ]
    
    try:
        subprocess.run(command, check=True)
        print(f"Canu command executed successfully for {fastq_file}.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running Canu for {fastq_file}: {e}")

def wtdbg2(fastq_file, base_output_dir, threads=24, genome_size="5m"):
    """
    执行 wtdbg2 命令，使用指定的 fastq 文件作为输入，并指定线程数和基因组大小
    在指定的 base_output_dir 下创建以 wtdbg2 为基础的输出文件夹
    """
    # 从文件路径中提取文件名（不包括扩展名）
    filename = os.path.basename(fastq_file)
    file_prefix = os.path.splitext(filename)[0]
    
    # 创建基础输出目录 'easynanometa_out/wtdbg2_out'
    wtdbg2_out_dir = os.path.join(base_output_dir, "wtdbg2_out")
    if not os.path.exists(wtdbg2_out_dir):
        os.makedirs(wtdbg2_out_dir)
        print(f"Created directory: {wtdbg2_out_dir}")
    
    # 在输出目录下创建以 file_prefix_wtdbg2 为前缀的子目录
    output_dir = os.path.join(wtdbg2_out_dir, f"{file_prefix}_wtdbg2", file_prefix)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")
    else:
        print(f"Directory already exists: {output_dir}")

    # 构建 wtdbg2 命令
    command = [
        "singularity", "exec", "easynanometa1_6.sif",
        "/bin/bash", "-c",
        f"source /tools/Miniconda3/bin/activate base && "
        f"/tools/Software/wtdbg2/wtdbg2.pl -t {threads} -x ont -g {genome_size} -o \"{output_dir}\" \"{fastq_file}\""
    ]
    
    # 打印命令行以供调试
    print("Running command:", " ".join(command))
    
    try:
        # 执行命令
        subprocess.run(command, check=True)
        print(f"wtdbg2 command executed successfully for {fastq_file}.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running wtdbg2 for {fastq_file}: {e}")

def nextdenove(fastq_file, base_output_dir, genome_size):
    """
    生成 {file_prefix}.fofn 文件，并为每个 .fofn 文件创建对应的 run.cfg 文件
    执行 NextDenovo 组装命令
    """
    filename = os.path.basename(fastq_file)
    file_prefix = os.path.splitext(filename)[0]
    
    nextdenove_out_dir = os.path.join(base_output_dir, "nextdenove_out")
    if not os.path.exists(nextdenove_out_dir):
        os.makedirs(nextdenove_out_dir)
        print(f"Created directory: {nextdenove_out_dir}")
    
    fofn_file = f"{file_prefix}.fofn"
    fofn_path = os.path.join(nextdenove_out_dir, fofn_file)
    
    # 生成 .fofn 文件
    with open(fofn_path, 'w') as f:
        f.write(f"{fastq_file}\n")
    
    # 为每个 .fofn 文件创建对应的 run.cfg 文件
    run_cfg_path = os.path.join(nextdenove_out_dir, f"{file_prefix}.cfg")
    with open(run_cfg_path, 'w') as f:
        f.write(f"""
###
[General]
job_type = local
job_prefix = nextDenovo
task = all # 'all', 'correct', 'assemble'
rewrite = yes # yes/no
deltmp = yes
rerun = 3
parallel_jobs = 2
input_type = raw
read_type = clr
input_fofn = {fofn_path}
workdir = ./{file_prefix}

[correct_option]
read_cutoff = 300bp #设置序列过滤长度
genome_size = {genome_size} #设置基因组组装大小
pa_correction = 2
sort_options = -m 1g -t 2
minimap2_options_raw =  -t 8
correction_options = -p 15

[assemble_option]
minimap2_options_cns =  -t 8
nextgraph_options = -a 1
###""")
    
    print(f"run.cfg file created at: {run_cfg_path}")
    
    # 执行 NextDenovo 组装命令
    command = (
        f"singularity exec easynanometa1_6.sif /bin/bash -c "
        f"\"source /tools/Miniconda3/bin/activate base && "
        f"/tools/Software/NextDenovo/nextDenovo {run_cfg_path}\""
    )
    
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"NextDenovo command executed successfully for {fastq_file}.")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while running NextDenovo for {fastq_file}: {e}")

def opera_ms(long_read_folder, short_read_folder, num_processors):
    """
    运行 OPERA-MS 软件的命令，使用 Singularity 容器来执行，并处理长读数据文件夹中的每个文件。

    参数:
    - long_read_folder: 长读数据文件夹的绝对路径
    - short_read_folder: 短读数据文件夹的绝对路径
    - num_processors: 使用的处理器数量
    """
    # 获取长读数据文件夹中的所有 .fastq 或 .fq 文件
    fastq_files = get_long_fastq_filenames(long_read_folder)
    
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
                f"--no-ref-clustering --no-polishing "
                f"--num-processors {num_processors} "
                f"--out-dir {output_dir}\""
            )

            # 执行命令
            subprocess.run(command, shell=True, check=True)
            print(f"OPERA-MS command executed with long-read file: {long_read_file} and output directory: {output_dir}")
    else:
        print("No fastq/fq files found in the long read folder.")
        

def spades_run(long_read_folder, short_read_folder, threads):
    """
    运行 SPAdes 的函数。根据指定的长读和短读数据文件夹执行命令。
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
            f'--meta -t {threads} -o {output_dir}"'
        )

        # 执行命令
        print(f"Running command: {command}")
        subprocess.run(command, shell=True)
def unicycler_run(long_read_folder, short_read_folder, threads):
    """
    运行 unicycler 的函数。根据指定的长读和短读数据文件夹执行命令。
    """
    # 获取所有长读数据文件
    long_read_files = get_long_fastq_filenames(long_read_folder)
    
    if not long_read_files:
        print("No .fastq or .fq files found in the long read folder.")
        return
    
    # 基础输出目录
    unicycler_out_dir = os.path.join(os.getcwd(), "easynanometa_out", "unicycler_out")
    
    if not os.path.exists(unicycler_out_dir):
        os.makedirs(unicycler_out_dir)
        print(f"Created base output directory: {unicycler_out_dir}")
    else:
        print(f"Base output directory already exists: {unicycler_out_dir}")
    
    # 依次处理每个长读文件
    for long_read_file in long_read_files:
        # 获取前缀名（假设前缀名为文件名的前部分，去掉扩展名）
        prefix = os.path.splitext(long_read_file)[0]
        
        # 为每个前缀创建单独的输出文件夹
        output_dir = os.path.join(unicycler_out_dir, prefix)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            print(f"Created output directory for prefix '{prefix}': {output_dir}")
        
        # 构建命令
        command = (
            f'singularity exec easynanometa1_6.sif /bin/bash -c "'
            f'source /tools/Miniconda3/bin/activate unicycler && '
            f'unicycler '
            f'-1 {short_read_folder}/{prefix}*1* '
            f'-2 {short_read_folder}/{prefix}*2* '
            f'-l {long_read_folder}/{long_read_file} '
            f'-t {threads} -o {output_dir}"'
        )

        # 执行命令
        print(f"Running command: {command}")
        subprocess.run(command, shell=True)        
        
def main():
    parser = argparse.ArgumentParser(description="Run Flye, Canu, Wtdbg2, NextDenove, OPERA-MS, SPAdes or Unicycler on fastq or fq files")
    subparsers = parser.add_subparsers(dest='tool', help='Tool to run: flye, canu, wtdbg2, nextdenove, OPERA-MS, SPAdes or Unicycler.')

    # Flye 子解析器
    flye_parser = subparsers.add_parser('flye', help='Run Flye tool')
    flye_parser.add_argument("-f", "--folder", required=True, help="Path to the folder containing fastq/fq files (absolute path)")
    flye_parser.add_argument("-t", "--threads", type=int, default=24, help="Number of threads (default: 24)")

    # Canu 子解析器
    canu_parser = subparsers.add_parser('canu', help='Run Canu tool')
    canu_parser.add_argument("-f", "--folder", required=True, help="Path to the folder containing fastq/fq files (absolute path)")
    canu_parser.add_argument("-g", "--genome-size", default="150m", help="Genome size (default: 150m)")

    # wtdbg2 子解析器
    wtdbg2_parser = subparsers.add_parser('wtdbg2', help='Run wtdbg2 tool')
    wtdbg2_parser.add_argument("-f", "--folder", required=True, help="Path to the folder containing fastq/fq files")
    wtdbg2_parser.add_argument("-t", "--threads", type=int, default=24, help="Number of threads to use for wtdbg2 (default: 24)")
    wtdbg2_parser.add_argument("-g", "--genome-size", default="5m", help="Genome size for wtdbg2 (default: 5m)")

    # NextDenovo 子解析器
    nextdenove_parser = subparsers.add_parser('nextdenove', help='Run NextDenovo tool')
    nextdenove_parser.add_argument("-f", "--folder", required=True, help="Path to the folder containing fastq/fq files (absolute path)")
    nextdenove_parser.add_argument("-g", "--genome-size", default="150m", help="Genome size (default: 150m)")
    
    # opera 子解析器
    opera_parser = subparsers.add_parser('opera', help='Run OPERA-MS tool')
    opera_parser.add_argument('-lf', '--long-read-folder', required=True, help="Absolute path to the long read data folder.")
    opera_parser.add_argument('-sf', '--short-read-folder', required=True, help="Absolute path to the short read data folder.")
    opera_parser.add_argument('-t', '--num-processors', type=int, default=24, help="Number of processors to use (default: 24).")
    

    # spades 子解析器
    spades_parser = subparsers.add_parser('spades', help='Run spades tool')
    spades_parser.add_argument('-lf', '--long_read_folder', type=str, required=True, help="Absolute path to the long read data folder.")
    spades_parser.add_argument('-sf', '--short_read_folder', type=str, required=True, help="Absolute path to the short read data folder.")
    spades_parser.add_argument('-t', '--threads', type=int, default=24, help="Number of processors to use (default: 24).")
    
    # unicycler 子解析器
    unicycler_parser = subparsers.add_parser('unicycler', help='Run unicycler tool')
    unicycler_parser.add_argument('-lf', '--long_read_folder', type=str, required=True, help="Absolute path to the long read data folder.")
    unicycler_parser.add_argument('-sf', '--short_read_folder', type=str, required=True, help="Absolute path to the short read data folder.")
    unicycler_parser.add_argument('-t', '--threads', type=int, default=48, help="Number of processors to use (default: 48).")
    
    args = parser.parse_args()
    
    if args.tool is None:
        parser.print_help()
        sys.exit(1)

    base_output_dir = create_base_output_directory()

    if args.tool == 'flye':
        filenames = get_fastq_filenames(args.folder)
        for fastq_file in filenames:
            flye(fastq_file, args.threads, base_output_dir)

    elif args.tool == 'canu':
        filenames = get_fastq_filenames(args.folder)
        for fastq_file in filenames:
            canu(fastq_file, base_output_dir, args.genome_size)

    elif args.tool == 'wtdbg2':
        filenames = get_fastq_filenames(args.folder)
        for fastq_file in filenames:
            wtdbg2(fastq_file, base_output_dir, args.threads, args.genome_size)

    elif args.tool == 'nextdenove':
        filenames = get_fastq_filenames(args.folder)
        for fastq_file in filenames:
            nextdenove(fastq_file, base_output_dir, args.genome_size)
    
    elif args.tool == 'opera':
        long_read_folder = args.long_read_folder
        short_read_folder = args.short_read_folder
        num_processors = args.num_processors
        fastq_files = get_long_fastq_filenames(long_read_folder)
        if fastq_files:
            for fastq_file in fastq_files:
                opera_ms(long_read_folder, short_read_folder, num_processors)       
    
    elif args.tool == 'spades':
        long_read_folder = args.long_read_folder
        short_read_folder = args.short_read_folder
        threads = args.threads
        fastq_files = get_long_fastq_filenames(long_read_folder)
        if fastq_files:
            for fastq_file in fastq_files:
                spades_run(long_read_folder, short_read_folder, threads) 
    elif args.tool == 'unicycler':
        long_read_folder = args.long_read_folder
        short_read_folder = args.short_read_folder
        threads = args.threads
        fastq_files = get_long_fastq_filenames(long_read_folder)
        if fastq_files:
            for fastq_file in fastq_files:
                unicycler_run(long_read_folder, short_read_folder, threads)
                
    else:
        print("No valid tool selected. Use -h for help.")
        sys.exit(1)

if __name__ == "__main__":
    main()
