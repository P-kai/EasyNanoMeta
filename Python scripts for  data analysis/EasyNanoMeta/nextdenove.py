import os
import subprocess
import argparse

def generate_fofn(file_name="long.fasta"):
    """生成 input.fofn 文件并将指定文件写入当前路径"""
    output_file = "input.fofn"  # 生成的文件在当前路径
    command = f"ls {file_name} > {output_file}"
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"{output_file} has been created with {file_name}.")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while generating {output_file}: {e}")

def get_fofn(file_name="input.fofn"):
    """获取 input.fofn 文件的绝对路径"""
    abs_path = os.path.abspath(file_name)
    if os.path.isfile(file_name):
        print(f"Absolute path of {file_name}: {abs_path}")
        return abs_path
    else:
        print(f"File {file_name} not found.")
        return None

def create_run_cfg(input_fofn_path, output_file="run.cfg"):
    """创建并生成 run.cfg 文件，使用绝对路径的 input.fofn 和固定的 workdir"""
    config_content = f"""
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
input_fofn = {input_fofn_path}
workdir = ./easynanometa_out/nextdenove_out

[correct_option]
read_cutoff = 300bp #设置序列过滤长度
genome_size = 308161 #设置基因组组装大小
pa_correction = 2
sort_options = -m 1g -t 2
minimap2_options_raw =  -t 8
correction_options = -p 15

[assemble_option]
minimap2_options_cns =  -t 8
nextgraph_options = -a 1
###
"""

    try:
        with open(output_file, 'w') as f:
            f.write(config_content.strip())
        print(f"{output_file} has been created successfully.")
        return os.path.abspath(output_file)
    except Exception as e:
        print(f"An error occurred while creating {output_file}: {e}")
        return None

def nextDenovo(run_cfg_path):
    """使用 Singularity 执行 NextDenovo 命令"""
    command = (f"singularity exec easynanometa1_6.sif /bin/bash -c "
               f"\"source /tools/Miniconda3/bin/activate base && /tools/Software/NextDenovo/nextDenovo {run_cfg_path}\"")
    try:
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(f"NextDenovo has been executed with {run_cfg_path}.")
        print("Standard Output:\n", result.stdout.decode())
        print("Standard Error:\n", result.stderr.decode())
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while executing NextDenovo: {e}")
        print("Standard Output:\n", e.stdout.decode())
        print("Standard Error:\n", e.stderr.decode())

def main():
    parser = argparse.ArgumentParser(description="Generate input.fofn, get its path, create run.cfg file, and execute NextDenovo.")
    parser.add_argument("-l", "--long_fasta", type=str, required=True, help="The file path of the long.fasta file")
    args = parser.parse_args()

    # Step 1: Generate input.fofn
    generate_fofn(args.long_fasta)

    # Step 2: Get absolute path of input.fofn
    input_fofn_abs_path = get_fofn("input.fofn")

    # Step 3: Create run.cfg using the absolute path of the input.fofn and fixed workdir
    if input_fofn_abs_path:
        run_cfg_path = create_run_cfg(input_fofn_abs_path)

        # Step 4: Execute NextDenovo using the generated run.cfg
        if run_cfg_path:
            nextDenovo(run_cfg_path)

if __name__ == "__main__":
    main()
