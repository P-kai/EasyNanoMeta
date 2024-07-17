import os
import argparse
import subprocess

def metaflye_run(flye_executable, raw_long_file, metaflye_threads):
    """
    运行 MetaFlye 的组装流程。
    """
    # 设置变量
    coa = f"{flye_executable} --meta"
    cob = "--nano-raw"
    coc = f"{raw_long_file}"
    cod = "--out-dir"
    coe = f"--threads {metaflye_threads}"

    # 获取文件名并去除后缀
    raw_long_name = os.path.splitext(os.path.basename(raw_long_file))[0]
    
    # 建立文件夹，名为 metaflye_out
    output_dir = 'metaflye_out'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")
    else:
        print(f"Directory already exists: {output_dir}")
   
    # 生成 code.sh 文件    
    with open('code.sh', 'w') as code_file:
            command = f'{coa} {cob} {coc} {cod} {output_dir}/{raw_long_name} {coe}\n'
            code_file.write(command)
    
    # 生成 activate.sh 文件
    with open('activate.sh', 'w') as activate_file:
        activate_file.write('#!/bin/bash\nchmod +x code.sh\n./code.sh\n')
    
    # 赋予 activate.sh 执行权限
    os.chmod('activate.sh', 0o755)
    
    # 执行 activate.sh 文件
    os.system('./activate.sh')
    
    # 删除 code.sh 文件
    os.remove('code.sh')
    
    print("MetaFlye assembly completed.")

    
def mv(raw_long_file):
    """
    修改metaflye_run结果文件名。
    """
    # 获取文件名并去除后缀
    raw_long_name = os.path.splitext(os.path.basename(raw_long_file))[0]  
    
    src = f'./metaflye_out/{raw_long_name}/assembly.fasta'
    dst = f'./metaflye_out/{raw_long_name}/{raw_long_name}.fasta'
    command = f'mv {src} {dst}'
    subprocess.run(command, shell=True, check=True)
    print(f"Moved {src} to {dst}.")

def sgs_lgs_fofn(raw_long_file, raw_short1_file, raw_short2_file):
    """
    生成 sgs.fofn 和 lgs.fofn 文件。
    """
    # 生成 sgs.fofn 文件
    sgs_output_file = 'sgs.fofn'
    with open(sgs_output_file, 'w') as sgs_file:
        sgs_file.write(f"{raw_short1_file}\n{raw_short2_file}\n")
    
    # 生成 lgs.fofn 文件
    lgs_output_file = 'lgs.fofn'
    with open(lgs_output_file, 'w') as lgs_file:
        lgs_file.write(f"{raw_long_file}\n")
    
    print(f"{sgs_output_file} and {lgs_output_file} generated.")

def create_run_cfg(raw_long_file):
    """
    创建并编辑 run.cfg 文件。
    """
    # 获取文件名并去除后缀
    raw_long_name = os.path.splitext(os.path.basename(raw_long_file))[0]
    
    cfg_content = f"""
[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 6
multithread_jobs = 5
genome = ./metaflye_out/{raw_long_name}/{raw_long_name}.fasta
genome_size = auto
workdir = ./short-and-long-reads-polish
polish_options = -p 8

[sgs_option]
sgs_fofn = ./sgs.fofn
sgs_options = -max_depth 100 -bwa

[lgs_option]
lgs_fofn = ./lgs.fofn
lgs_options = -min_read_len 1k -max_depth 100
lgs_minimap2_options = -x map-ont
    """
    with open('run.cfg', 'w') as cfg_file:
        cfg_file.write(cfg_content)
    print("run.cfg created.")

def nextpolish(executable):
    """
    执行 NextPolish 命令的函数。
    """
    command = f'{executable} run.cfg'
    subprocess.run(command, shell=True, check=True)
    print("NextPolish executed successfully.")

def mv1(raw_long_file):
    """
    修改nextpolish结果文件名。
    """
    # 获取文件名并去除后缀
    raw_long_name = os.path.splitext(os.path.basename(raw_long_file))[0]
    
    src = f'./short-and-long-reads-polish/genome.nextpolish.fasta'
    dst = f'./short-and-long-reads-polish/{raw_long_name}.nextpolish.fasta'
    command = f'mv {src} {dst}'
    subprocess.run(command, shell=True, check=True)
    print(f"Moved {src} to {dst}.")
    
def minimap2(raw_long_file):
    """
    对组装基因组进行索引创建
    """
    # 获取文件名并去除后缀
    raw_long_name = os.path.splitext(os.path.basename(raw_long_file))[0]
    
    command = f"minimap2 -d catalogue.mmi ./short-and-long-reads-polish/{raw_long_name}.nextpolish.fasta"
    subprocess.run(command, shell=True, check=True)
   

def bam(raw_long_file):
    """
    比对获取bam文件
    """
    # 创建bam文件夹
    os.makedirs("bam", exist_ok=True)
    print("The bam folder was created successfully.")
    
    # 获取文件名并去除后缀
    raw_long_name = os.path.splitext(os.path.basename(raw_long_file))[0]
    
    command = f"minimap2 -t 8 -N 5 -ax map-ont catalogue.mmi {raw_long_file} | samtools view -F 3584 -b --threads 8 > ./bam/{raw_long_name}.bam"
    subprocess.run(command, shell=True, check=True)
    print("BAM file generated.")
    

def samtools(raw_long_file):
    """
    对bam文件进行排序
    """
    # 获取文件名并去除后缀
    raw_long_name = os.path.splitext(os.path.basename(raw_long_file))[0]
    
    command = f"samtools sort -@ 10 ./bam/{raw_long_name}.bam > ./bam/{raw_long_name}.sorted.bam"
    subprocess.run(command, shell=True, check=True)
    print("BAM file sorted.")

def semi_bin(semi_bin_path, raw_long_file):
    """
    使用SemiBin进行基因组分箱
    """
    # 获取文件名并去除后缀
    raw_long_name = os.path.splitext(os.path.basename(raw_long_file))[0]
    
    activate_conda_command = "source /home/tools_pk/miniconda3/etc/profile.d/conda.sh && conda activate /home/tools_pk/miniconda3/envs/SemiBin"
    semi_bin_command = f"{semi_bin_path} single_easy_bin -i ./short-and-long-reads-polish/{raw_long_name}.nextpolish.fasta --sequencing-type long_read -b ./bam/{raw_long_name}.sorted.bam -o sembin_out --environment global"
    deactivate_conda_command = "conda deactivate"
    
    command = f"{activate_conda_command} && {semi_bin_command} && {deactivate_conda_command}"
    subprocess.run(command, shell=True, check=True, executable='/bin/bash')
    print("SemiBin command executed.")

def checkm2(checkm2_path, checkm2_database):
    """
    使用CheckM2预测宏基因组组装基因组的完整性和污染
    """
    activate_conda_command = "source /home/tools_pk/miniconda3/etc/profile.d/conda.sh && conda activate /home/tools_pk/miniconda3/envs/checkm2"
    checkm2_command = f"{checkm2_path} predict --database_path {checkm2_database} --threads 30 --input ./sembin_out/output_bins/* --output-directory ./checkm2_out"
    deactivate_conda_command = "conda deactivate"
    
    command = f"{activate_conda_command} && {checkm2_command} && {deactivate_conda_command}"
    try:
        subprocess.run(command, shell=True, check=True, executable='/bin/bash')
        print("The checkm2 command ran successfully.")
    except subprocess.CalledProcessError as e:
        print(f"CheckM2 command failed to run: {e}")
    except Exception as e:
        print(f"An error occurred while running CheckM2: {e}")

def gtdbtk(gtdbtk_path, gtdbtk_database):
    """
    使用GTDB-Tk对基因组进行分类
    """
    activate_conda_command = "source /home/tools_pk/miniconda3/etc/profile.d/conda.sh && conda activate /home/tools_pk/miniconda3/envs/gtdbtk"
    export_command = f"export GTDBTK_DATA_PATH={gtdbtk_database}"
    gtdbtk_command = f"{gtdbtk_path} classify_wf --genome_dir ./sembin_out/output_bins/ --extension fa --skip_ani_screen --out_dir ./gtdbtk"
    deactivate_conda_command = "conda deactivate"
    
    command = f"{activate_conda_command} && {export_command} && {gtdbtk_command} && {deactivate_conda_command}"
    try:
        subprocess.run(command, shell=True, check=True, executable='/bin/bash')
        print("The gtdbtk command ran successfully.")
    except subprocess.CalledProcessError as e:
        print(f"GTDB-Tk command failed to run: {e}")
    except Exception as e:
        print(f"An error occurred while running GTDB-Tk: {e}")

def main():
    """
    主程序入口。
    """
    parser = argparse.ArgumentParser(description="Main script to run multiple bioinformatics tools.")
    parser.add_argument('--flye', default='/Tools/software/Flye-2021-2.9/bin/flye', 
                        help='Path to Flye executable (default: /Tools/software/Flye-2021-2.9/bin/flye).')
    parser.add_argument('--threads', type=int, default=24, help='Number of threads for MetaFlye.')   
    parser.add_argument('--rawlong', required=True, help='Input raw long-read file.')
    parser.add_argument('--rawshort1', required=True, help='Input raw short-read1 file.')
    parser.add_argument('--rawshort2', required=True, help='Input raw short-read2 file.')
    parser.add_argument('--nextpolish', default='/home/tools_pk/tools/NextPolish/nextPolish', 
                        help='Path to NextPolish executable (default: /home/tools_pk/tools/NextPolish/nextPolish)')
    parser.add_argument('--semibin', default='/home/tools_pk/miniconda3/envs/SemiBin/bin/SemiBin', 
                        help='Path to SemiBin executable (default: /home/tools_pk/miniconda3/envs/SemiBin/bin/SemiBin).')
    parser.add_argument("--checkm2", type=str, default="/home/tools_pk/miniconda3/envs/checkm2/bin/checkm2",
                        help="Path to the CheckM2 executable (default: /home/tools_pk/miniconda3/envs/checkm2/bin/checkm2).")
    parser.add_argument("--checkm2_database", type=str, default="/home/tools_pk/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd",
                        help="Path to the CheckM2 database (default: /home/tools_pk/databases/checkm2/CheckM2_database/uniref100.KO.1.dmnd).")
    parser.add_argument("--gtdbtk", type=str, default="/home/tools_pk/miniconda3/envs/gtdbtk/bin/gtdbtk",
                        help="Path to the GTDB-Tk executable (default: /home/tools_pk/miniconda3/envs/gtdbtk/bin/gtdbtk).")
    parser.add_argument("--gtdbtk_database", type=str, default="/backup/database/gtdbtk/release214",
                        help="Path to the GTDB-Tk database (default: /backup/database/gtdbtk/release214).")

    args = parser.parse_args()

    metaflye_run(args.flye, args.rawlong, args.threads)
    mv(args.rawlong)
    create_run_cfg(args.rawlong)
    sgs_lgs_fofn(args.rawlong, args.rawshort1, args.rawshort2)
    nextpolish(args.nextpolish)
    mv1(args.rawlong)
    minimap2(args.rawlong)
    bam(args.rawlong)
    samtools(args.rawlong)
    semi_bin(args.semibin, args.rawlong)
    checkm2(args.checkm2, args.checkm2_database)
    gtdbtk(args.gtdbtk, args.gtdbtk_database)

if __name__ == "__main__":
    main()
