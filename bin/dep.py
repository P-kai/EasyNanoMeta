import os
import subprocess
import sys
import shutil 
import pandas as pd

# Flye functionality
def run_flye(input_dir, output_dir, singularity_image, flye_executable, threads=24):
    # Check if input directory exists
    if not os.path.exists(input_dir):
        raise ValueError(f"Input directory {input_dir} does not exist")
    
    # Check if output directory exists, if not, create it
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Iterate through all fastq files in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.fastq'):
            input_file = os.path.join(input_dir, file_name)
            sample_name = os.path.splitext(file_name)[0]  # Get file name without extension
            output_path = os.path.join(output_dir, f"{sample_name}_flye")

            # Create output path if it does not exist
            if not os.path.exists(output_path):
                os.makedirs(output_path)

            # Generate Singularity Flye command
            command = [
                "singularity", "exec", "-B",
                f"{input_dir}:/mnt",
                singularity_image,
                flye_executable,
                "--nano-raw", f"/mnt/{file_name}",
                "--threads", str(threads),
                "--meta",
                "-o", output_path
            ]

            # Print the command to verify it
            print(f"Running Flye on {file_name} with output to {output_path}")
            print("Command:", " ".join(command))

            # Execute the command
            try:
                subprocess.run(command, shell=True, check=True)
                print(f"Flye assembly for {file_name} completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error running Flye on {file_name}: {e}")

# Kraken2 functionality
def run_kraken2(input_dir, output_dir, singularity_image, kraken2_executable, kraken2_db, threads=24):
    # Check if input directory exists
    if not os.path.exists(input_dir):
        raise ValueError(f"Input directory {input_dir} does not exist")
    
    # Check if output directory exists, if not, create it
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Iterate through all fastq files in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith('.fastq'):
            input_file = os.path.join(input_dir, file_name)
            sample_name = os.path.splitext(file_name)[0]  # Get file name without extension

            # Generate output file paths
            report_file = os.path.join(output_dir, f"{sample_name}_kraken2_report")
            result_file = os.path.join(output_dir, f"{sample_name}_kraken2_result")

            # Generate Singularity Kraken2 command
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

            # Print the command to verify it
            print(f"Running Kraken2 on {file_name} with report output to {report_file}")
            print("Command:", " ".join(command))

            # Execute the command
            try:
                subprocess.run(command, check=True)
                print(f"Kraken2 analysis for {file_name} completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error running Kraken2 on {file_name}: {e}")


# Abricate functionality
def run_abricate(input_dir, output_dir, singularity_image, abricate_executable, db='ncbi'):
    if not os.path.exists(input_dir):
        raise ValueError(f"Input directory {input_dir} does not exist")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for file_name in os.listdir(input_dir):
        if file_name.endswith('.fa'):
            input_file = os.path.join(input_dir, file_name)
            sample_name = os.path.splitext(file_name)[0]
            output_file = os.path.join(output_dir, f"{sample_name}_abricate_output.txt")

            command = [
                "singularity", "exec", "-B",
                f"{input_dir}:/mnt",
                singularity_image,
                "/bin/bash", "-c",
                f'"source /tools/Miniconda3/bin/activate abricate && '
                f'{abricate_executable} --db {db} /mnt/{file_name} > {output_file}"'
            ]
            print(f"Running Abricate on {file_name} with output to {output_file}")
            try:
                subprocess.run(" ".join(command), shell=True, check=True)
                print(f"Abricate analysis for {file_name} completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error running Abricate on {file_name}: {e}")

def run_adapters_removal(input_dir, output_dir, singularity_image, adapters_removal_executable, num_threads):
    if not os.path.exists(input_dir):
        raise ValueError(f"Input directory {input_dir} does not exist")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for file_name in os.listdir(input_dir):
        if file_name.endswith('.fastq'):
            input_file = os.path.join(input_dir, file_name)
            sample_name = os.path.splitext(file_name)[0]

            command = (
                f"singularity exec -B {input_dir}:/mnt {singularity_image} /bin/bash -c \"source /tools/Miniconda3/bin/activate porechop_abi && "
                f"{adapters_removal_executable} --ab_initio -i /mnt/{file_name} "
                f"-o {output_dir}/{sample_name}_output.fastq "
                f"-t {num_threads} \""
            )

            print(f"Running adapters removal for {file_name} with output to {output_dir}/{sample_name}_output.fastq")
            print("Command:", " ".join(command))

            try:
                subprocess.run(command, shell=True, check=True)
                print(f"Adapters removal for {file_name} completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error occurred while processing {file_name}: {e}")

def run_host_removal(host_removal_reference, input_dir, output_dir, singularity_image, num_threads):
    if not os.path.exists(input_dir):
        raise ValueError(f"Input directory {input_dir} does not exist")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    reference_prefix = os.path.splitext(os.path.basename(host_removal_reference))[0]

    fasta_files = []  

    for file_name in os.listdir(input_dir):
        if file_name.endswith("_output.fastq"):
            input_file = os.path.join(input_dir, file_name)
            sample_name = os.path.splitext(file_name)[0]
            fasta_output_path = os.path.join(output_dir, f"{sample_name.replace('_output', '')}.fasta")

            command = (
                f"singularity exec {singularity_image} /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"seqkit fq2fa {input_file} > {fasta_output_path}\""
            )
            
            try:
                subprocess.run(command, shell=True, check=True)
                print(f"Converted fastq to fasta for {file_name} completed successfully.")
                fasta_files.append(fasta_output_path)
            except subprocess.CalledProcessError as e:
                print(f"Error occurred while processing {file_name}: {e}")

    hostremoval_reference_directory = os.path.dirname(host_removal_reference)
    fasta_file_name = os.path.basename(host_removal_reference)
    
    minimap2_index_cmd = f"singularity exec -B {hostremoval_reference_directory}:/mnt {singularity_image} /bin/bash -c \"source /tools/Miniconda3/bin/activate host_removal && minimap2 -d {output_dir}/{reference_prefix}.min /mnt/{fasta_file_name}\""
 
    try:
        print(f"Building minimap2 index for host genome {host_removal_reference}...")
        subprocess.run(minimap2_index_cmd, shell=True, check=True)
        print(f"Minimap2 index created: {output_dir}/{reference_prefix}.min")
    except subprocess.CalledProcessError as e:
        print(f"Error during minimap2 index creation: {e}")

    for fasta_file in fasta_files:
        file_prefix = os.path.splitext(os.path.basename(fasta_file))[0]

        minimap2_align_cmd = f"singularity exec {singularity_image} /bin/bash -c \"source /tools/Miniconda3/bin/activate host_removal && minimap2 -ax map-ont -t {num_threads} {output_dir}/{reference_prefix}.min {fasta_file} -o {output_dir}/{file_prefix}_minimap.sam\""
        try:
            print(f"Aligning {fasta_file} to host genome...")
            subprocess.run(minimap2_align_cmd, shell=True, check=True)
            print(f"Alignment completed: {output_dir}/{file_prefix}_minimap.sam")
        except subprocess.CalledProcessError as e:
            print(f"Error during alignment: {e}")
            continue

        extract_unmaped_cmd = f"singularity exec -B {hostremoval_reference_directory}:/mnt {singularity_image} /bin/bash -c \"source /tools/Miniconda3/bin/activate host_removal && samtools view -bS -@ {num_threads} -T /mnt/{fasta_file_name} -f 4 {output_dir}/{file_prefix}_minimap.sam > {output_dir}/{file_prefix}_unmaped_minimap.bam\""
        try:
            print(f"Extracting unmaped reads from {output_dir}/{file_prefix}_minimap.sam...")
            subprocess.run(extract_unmaped_cmd, shell=True, check=True)
            print(f"Unmaped reads extracted: {output_dir}/{file_prefix}_unmaped_minimap.bam")
        except subprocess.CalledProcessError as e:
            print(f"Error during unmaped reads extraction: {e}")
            continue

        sort_bam_cmd = f"singularity exec {singularity_image} /bin/bash -c \"source /tools/Miniconda3/bin/activate host_removal && samtools sort -n {output_dir}/{file_prefix}_unmaped_minimap.bam -o {output_dir}/{file_prefix}_unmaped_sorted_minimap.bam\""
        try:
            print(f"Sorting BAM files...")
            subprocess.run(sort_bam_cmd, shell=True, check=True)
            print(f"Sorted BAM files: {output_dir}/{file_prefix}_unmaped_sorted_minimap.bam")
        except subprocess.CalledProcessError as e:
            print(f"Error during BAM file sorting: {e}")
            continue
        
        bam_to_fastq_cmd = f"singularity exec {singularity_image} /bin/bash -c \"source /tools/Miniconda3/bin/activate host_removal && bedtools bamtofastq -i {output_dir}/{file_prefix}_unmaped_sorted_minimap.bam -fq {output_dir}/{file_prefix}_fitted_raw.fastq\""
        try:
            print(f"Converting BAM files to FASTQ files...")
            subprocess.run(bam_to_fastq_cmd, shell=True, check=True)
            print(f"Conversion completed: {output_dir}/{file_prefix}_fitted_raw.fastq")
        except subprocess.CalledProcessError as e:
            print(f"Error during BAM to FASTQ conversion: {e}")
            continue

def run_centrifuge(centrifuge_executable, kreport_executable, centrifuge_db, input_dir, output_dir, singularity_image, num_threads):
    if not os.path.exists(input_dir):
        raise ValueError(f"Input directory {input_dir} does not exist")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    for file_name in os.listdir(input_dir):
        if file_name.endswith("_fitted_raw.fastq"):
            sample_name = os.path.splitext(file_name)[0]
            input_file = os.path.join(input_dir, file_name)
            centrifuge_report_path = os.path.join(output_dir, f"{sample_name.replace('_fitted_raw', '')}_report")
            centrifuge_classification_path = os.path.join(output_dir, f"{sample_name.replace('_fitted_raw', '')}_result")
                   
            centrifuge_db_directory = os.path.dirname(centrifuge_db)
            db_name = os.path.basename(centrifuge_db)
            centrifuge_cmd = (
                f"singularity exec -B {centrifuge_db_directory}:/mnt {singularity_image} {centrifuge_executable} "
                f"-p {num_threads} -x /mnt/{db_name} "
                f"-q {input_file} "
                f"--report-file {centrifuge_report_path} "
                f"-S {centrifuge_classification_path}"
            )   
            try:
                print(f"Running Centrifuge for {file_name}...")
                subprocess.run(centrifuge_cmd, shell=True, check=True)
                print(f"Centrifuge classification completed for {file_name}")
            except subprocess.CalledProcessError as e:
                print(f"Error during Centrifuge classification for {file_name}: {e}")
                continue
    
            centrifuge_to_kraken_cmd = (
                f"singularity exec -B {centrifuge_db_directory}:/mnt {singularity_image} {kreport_executable} "
                f"-x /mnt/{db_name} "
                f"{centrifuge_classification_path} > "
                f"{output_dir}/{sample_name}_kraken_report"
            )
            
            try:
                print(f"Converting Centrifuge result to Kraken format for {file_name}...")
                subprocess.run(centrifuge_to_kraken_cmd, shell=True, check=True)
                print(f"Centrifuge result converted to Kraken report for {file_name}")
            except subprocess.CalledProcessError as e:
                print(f"Error during Centrifuge-to-Kraken conversion for {file_name}: {e}")
                continue

def run_arg_abundance(input_dir, output_dir, singularity_image, fasta_dir):
    def get_file_prefix(filename, suffixes):
        for suffix in suffixes:
            if filename.endswith(suffix):
                return filename[:-len(suffix)]
        return filename

    def bp_to_gb(bp):
        try:
            return float(bp) / 1e9
        except ValueError:
            print(f"Invalid genome size: {bp}")
            return None

    def get_genome_size_from_file(file_path):
        try:
            df = pd.read_csv(file_path, sep='\t', header=None)
            genome_size = df.iloc[1, 4] 
            return genome_size
        except Exception as e:
            print(f"Error reading genome size from {file_path}: {e}")
            return None
    
    if not os.path.exists(input_dir):
        raise ValueError(f"Input directory {input_dir} does not exist")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fasta_files = []
    for fasta_file in os.listdir(fasta_dir):
        if fasta_file.endswith(".fasta") or fasta_file.endswith(".fa"):
            fasta_path = os.path.join(fasta_dir, fasta_file)
            fasta_files.append(fasta_path)

    genome_sizes = {}
    for fasta_file in fasta_files:
        file_prefix = get_file_prefix(os.path.basename(fasta_file), ['.fasta', '.fa'])
        seqkit_cmd = (
            f"singularity exec {singularity_image} /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
            f"seqkit stats {fasta_file} --all --tabular > {output_dir}/{file_prefix}_genome.csv\""
        )
        try:
            subprocess.run(seqkit_cmd, shell=True, check=True)
            excel_file = os.path.join(output_dir, f"{file_prefix}_genome.csv")
            genome_size = get_genome_size_from_file(excel_file)
            if genome_size is not None:
                genome_sizes[file_prefix] = genome_size
                print(f"Genome size of {fasta_file}: {genome_size} bp")
            else:
                print(f"Could not find genome size for {fasta_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error getting genome size for {fasta_file}: {e}")
            continue

    for ncbi_result_file in os.listdir(input_dir):
        if ncbi_result_file.endswith("_ncbi_result"):
            file_prefix = get_file_prefix(ncbi_result_file, ['_ncbi_result'])
            ncbi_result_path = os.path.join(input_dir, ncbi_result_file)
        
            # 检查文件是否为空
            if os.path.getsize(ncbi_result_path) == 0:
                continue  # 如果文件为空，则跳过后续处理
        
            abundance_output_path = os.path.join(output_dir, f"{file_prefix}_abundance_result")
        
            genome_size_bp = genome_sizes.get(file_prefix, "unknown")
            genome_size_gb = bp_to_gb(genome_size_bp)
        
            if genome_size_gb is None:
                continue  # 跳过无效的基因组大小

            genome_size_gb_value = f"{genome_size_gb:.3f}"

            abundance_calculate_cmd = (
                f"singularity exec {singularity_image} /bin/bash -c \"source /tools/Miniconda3/bin/activate base && python "
                f"/tools/Software/abundance_calculate.py --i {ncbi_result_path} --data_size {genome_size_gb_value} -p {file_prefix} -o {output_dir}\""
            )
            print(f"Calculating abundance for {ncbi_result_file} with genome size {genome_size_gb_value} Gb ...")
            subprocess.run(abundance_calculate_cmd, shell=True, check=True)
            print(f"Abundance calculation completed: {abundance_output_path}")

def run_nextpolish(input_dir, output_dir, singularity_image, num_threads, flye_dir):
    if not os.path.exists(input_dir):
        raise ValueError(f"Input directory {input_dir} does not exist")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for fastq_file in os.listdir(input_dir):
        if fastq_file.endswith("_fitted_raw.fastq"):
            file_prefix = fastq_file.replace("_fitted_raw.fastq", "")
            fastq_path = os.path.join(input_dir, fastq_file)
            lgs_out_path = os.path.abspath(os.path.join(output_dir, f"{file_prefix}.fofn"))  # 使用绝对路径

            # 创建 lgs_out_path 文件
            try:
                with open(lgs_out_path, 'w') as lgs_file:
                    lgs_file.write(fastq_path + "\n")
                print(f"fofn file created successfully for {fastq_file}. Output file: {lgs_out_path}")
            except Exception as e:
                print(f"Error during fofn file creation for {fastq_file}: {e}")
                continue

            cfg_file_path = os.path.abspath(os.path.join(output_dir, f"{file_prefix}.cfg"))  # 使用绝对路径
            genome_path = os.path.abspath(os.path.join(flye_dir, f"{file_prefix}.fasta"))  # 使用绝对路径
            workdir_path = os.path.abspath(os.path.join(output_dir, f"{file_prefix}_nextpolish_out"))  # 使用绝对路径
            
            # 创建配置文件内容
            cfg_content = f"""\

[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 6
multithread_jobs = {num_threads}
genome = {genome_path} #组装结果文件
genome_size = auto
workdir = {workdir_path}
polish_options = -p {num_threads}

[lgs_option]
lgs_fofn = {lgs_out_path}
lgs_options = -min_read_len 1k -max_depth 100
lgs_minimap2_options = -x map-ont
"""

            # 创建 cfg 文件
            try:
                with open(cfg_file_path, 'w') as cfg_file:
                    cfg_file.write(cfg_content)
                print(f"cfg file created successfully. Output file: {cfg_file_path}")
            except Exception as e:
                print(f"Error during cfg file creation for {lgs_out_path}: {e}")
                continue  # 如果创建失败，继续下一个文件

            # 执行 NextPolish 命令
            nextpolish_cmd = (
                f"singularity exec {singularity_image} /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"/tools/Software/NextPolish/nextPolish {cfg_file_path}\""
            )
            
            try:
                print(f"Performing calibration procedures for {file_prefix}.fasta...")
                subprocess.run(nextpolish_cmd, shell=True, check=True)
                print(f"Calibration procedures completed successfully for {file_prefix}.fasta.")
            except subprocess.CalledProcessError as e:
                print(f"Error during calibration procedures for {file_prefix}.fasta: {e}")
                continue  # 如果出错，继续下一个文件

def run_semi_bin(input_dir, output_dir, singularity_image, num_threads, host_removal_dir):
    if not os.path.exists(input_dir):
        raise ValueError(f"Input directory {input_dir} does not exist")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    bam_out_folder = os.path.join(output_dir, "bam_out")
    
    if not os.path.exists(bam_out_folder):
        try:
            os.makedirs(bam_out_folder)
            print(f"Created folder: {bam_out_folder}")
        except OSError as e:
            print(f"Error creating folder {bam_out_folder}: {e}")
            sys.exit(1)
    else:
        print(f"Folder already exists: {bam_out_folder}")

    prefixes = []

    for folder_name in os.listdir(input_dir):
        folder_path = os.path.join(input_dir, folder_name)
        
        if os.path.isdir(folder_path) and folder_name.endswith("_nextpolish_out"):
            prefix = folder_name.rsplit("_nextpolish_out", 1)[0]
            prefixes.append((folder_path, prefix))
    
    for folder_path, prefix in prefixes:
        fasta_files = []
        
        for file_name in os.listdir(folder_path):
            if file_name == f"{prefix}_nextpolish.fasta":
                fasta_files.append(os.path.join(folder_path, file_name))
        
        for fasta_file in fasta_files:
            minimap2_index_cmd = (
                f"singularity exec {singularity_image} /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"minimap2 -d {os.path.join(folder_path, 'catalogue.mmi')} {fasta_file}\""
            )
            try:
                print(f"Creating index file for {fasta_file}...")
                subprocess.run(minimap2_index_cmd, shell=True, check=True)
                print(f"Index procedures completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error during index procedures for {fasta_file}: {e}")
                continue

    for fastq_file in os.listdir(host_removal_dir):
        if fastq_file.endswith("_fitted_raw.fastq"):
            file_prefix = fastq_file.replace("_fitted_raw.fastq", "")
            fastq_path = os.path.join(host_removal_dir, fastq_file)
            bam_out_path = os.path.join(bam_out_folder, f"{file_prefix}.bam")
            
            minimap2_cmd = (
                f"singularity exec {singularity_image} /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"minimap2 -t {num_threads} -N 5 -ax map-ont {os.path.join(folder_path, 'catalogue.mmi')} {fastq_path} | "
                f"samtools view -F 3584 -b --threads 8 > {bam_out_path}\""
            )
            try:
                print(f"BAM file generated for {fastq_path}...")
                subprocess.run(minimap2_cmd, shell=True, check=True)
                print(f"BAM file generated successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error during BAM file generation for {fastq_path}: {e}")
                continue
            
            sorted_bam_out_path = os.path.join(bam_out_folder, f"{file_prefix}.sorted.bam")
            samtools_cmd = (
                f"singularity exec {singularity_image} /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"samtools sort -@ 10 {bam_out_path} > {sorted_bam_out_path}\""
            )
            try:
                print(f"BAM file sorted for {file_prefix}.bam...")
                subprocess.run(samtools_cmd, shell=True, check=True)
                print(f"BAM file sorted successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error during BAM file sorting for {file_prefix}.bam: {e}")
                continue
    
    for folder_path, prefix in prefixes:
        fasta_files = []
        
        for file_name in os.listdir(folder_path):
            if file_name == f"{prefix}_nextpolish.fasta":
                fasta_files.append(os.path.join(folder_path, file_name))
    
        for fasta_file in fasta_files: 
            semi_bin_cmd = (
                f"singularity exec {singularity_image} /bin/bash -c \"source /tools/Miniconda3/bin/activate SemiBin && "
                f"SemiBin single_easy_bin -i {fasta_file} --sequencing-type long_read -b {sorted_bam_out_path} "
                f"-o {os.path.join(output_dir, f'{prefix}_bin_out')} --environment global\""
            )
            try:
                print(f"Performing bin for {fasta_file}...")
                subprocess.run(semi_bin_cmd, shell=True, check=True)
                print(f"Perform bin successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error during bin procedures for {fasta_file}: {e}")
                continue

def run_checkm2(input_dir, output_dir, singularity_image, num_threads, checkm2_db):
    if not os.path.exists(input_dir):
        raise ValueError(f"Input directory {input_dir} does not exist")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    prefixes = []
    
    for folder_name in os.listdir(input_dir):
        folder_path = os.path.join(input_dir, folder_name)
        
        if os.path.isdir(folder_path) and folder_name.endswith("_bin_out"):
            prefix = folder_name.rsplit("_bin_out", 1)[0]
            prefixes.append((folder_path, prefix))

            for subfolder_name in os.listdir(folder_path):
                output_bins_folder = os.path.join(folder_path, subfolder_name)
                
                if os.path.isdir(output_bins_folder) and subfolder_name == "output_bins":
                    print(f"Processing {output_bins_folder} for prefix {prefix}")

                    checkm2_db_directory = os.path.dirname(checkm2_db)
                    db_file_name = os.path.basename(checkm2_db)

                    checkm2_cmd = (
                        f"singularity exec -B {checkm2_db_directory}:/mnt {singularity_image} /bin/bash -c \"source /tools/Miniconda3/bin/activate checkm2 && "
                        f"checkm2 predict --database_path /mnt/{db_file_name} --threads {num_threads} --input {output_bins_folder}/* --output-directory {output_dir}/{prefix}_checkm2_out\""
                    )
                    try:
                        print(f"Perform checkm2 for {folder_path}...")
                        subprocess.run(checkm2_cmd, shell=True, check=True)
                        print(f"Perform checkm2 successfully.")
                    except subprocess.CalledProcessError as e:
                        print(f"Error during checkm2 procedures for {folder_path}: {e}")
                        continue

def run_gtdbtk(input_dir, output_dir, singularity_image, num_threads, gtdbtk_db):
    if not os.path.exists(input_dir):
        raise ValueError(f"Input directory {input_dir} does not exist")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    prefixes = []
    
    for folder_name in os.listdir(input_dir):
        folder_path = os.path.join(input_dir, folder_name)
        
        if os.path.isdir(folder_path) and folder_name.endswith("_bin_out"):
            prefix = folder_name.rsplit("_bin_out", 1)[0]
            prefixes.append((folder_path, prefix))

            for subfolder_name in os.listdir(folder_path):
                output_bins_folder = os.path.join(folder_path, subfolder_name)
                
                if os.path.isdir(output_bins_folder) and subfolder_name == "output_bins":
                    print(f"Processing {output_bins_folder} for prefix {prefix}")

                    gtdbtk_cmd = (
                        f"singularity exec -B {gtdbtk_db}:/mnt {singularity_image} /bin/bash -c \"source /tools/Miniconda3/bin/activate gtdbtk-2.2.6 && export GTDBTK_DATA_PATH=/mnt && "
                        f"gtdbtk classify_wf --genome_dir {output_bins_folder} --extension fa --skip_ani_screen --out_dir {output_dir}/{prefix}_gtdbtk_out\""
                    )
                    try:
                        print(f"Perform gtdbtk for {folder_path}...")
                        subprocess.run(gtdbtk_cmd, shell=True, check=True)
                        print(f"Perform gtdbtk successfully.")
                    except subprocess.CalledProcessError as e:
                        print(f"Error during gtdbtk procedures for {folder_path}: {e}")
                        continue
