#!/usr/bin/env python
################################################################
#Authors: Peng Kai; Li changan
#Email: 008719@yzu.edu.cn
#
#This is main script for EasyNanaMeta pipeline
################################################################
import os
import argparse
import subprocess
import sys
import shutil 
import pandas as pd

def create_output_folder():
    output_folder = os.path.join(os.getcwd(), "easynanometa_result")
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Created folder: {output_folder}")
    else:
        print(f"Folder already exists: {output_folder}")
    
    adapters_removal_out_folder = os.path.join(output_folder, "adapters_removal_out")
    if not os.path.exists(adapters_removal_out_folder):
        os.makedirs(adapters_removal_out_folder)
        print(f"Created folder: {adapters_removal_out_folder}")
    else:
        print(f"Folder already exists: {adapters_removal_out_folder}")
    
    return output_folder, adapters_removal_out_folder

def find_fastq_files(folder_path):
    fastq_files = {}
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith((".fastq", ".fq")):
                prefix = os.path.splitext(file)[0]
                fastq_files[prefix] = os.path.join(root, file)
    return fastq_files

def adapters_removal(file_path, file_prefix, output_folder, num_threads, sif_path):
    adapters_removal_out_folder = os.path.join(output_folder, "adapters_removal_out")
    
    command = f"singularity exec {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate porechop_abi && porechop_abi --ab_initio -i {file_path} -o {adapters_removal_out_folder}/{file_prefix}_output.fastq -t {num_threads}\""
    try:
        print(f"Running adapters removal for {file_path} with {num_threads} threads...")
        subprocess.run(command, shell=True, check=True)
        print(f"Adapters removal completed for {file_path}")
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while processing {file_prefix}: {e}")

def host_removal(host_removal_reference, num_threads, sif_path):
    host_removal_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "host_removal_out")
    if not os.path.exists(host_removal_out_folder):
        os.makedirs(host_removal_out_folder)
        print(f"Created folder: {host_removal_out_folder}")
    else:
        print(f"Folder already exists: {host_removal_out_folder}")

    reference_prefix = os.path.splitext(os.path.basename(host_removal_reference))[0]

    adapters_removal_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "adapters_removal_out")
    fasta_files = []  

    for fastq_file in os.listdir(adapters_removal_out_folder):
        if fastq_file.endswith("_output.fastq"):
            file_prefix = os.path.splitext(fastq_file)[0]  
            fastq_path = os.path.join(adapters_removal_out_folder, fastq_file)
            fasta_output_path = os.path.join(host_removal_out_folder, f"{file_prefix.replace('_output', '')}.fasta")
            fastq_to_fasta_cmd = (
                f"singularity exec {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"seqkit fq2fa {fastq_path} > {fasta_output_path}\""
            )
            try:
                subprocess.run(fastq_to_fasta_cmd, shell=True, check=True)
                print(f"Converted {fastq_path} to {fasta_output_path}")
                fasta_files.append(fasta_output_path)
            except subprocess.CalledProcessError as e:
               print(f"Error converting {fastq_path} to {fasta_output_path}: {e}")  
               
    hostremoval_reference_directory = os.path.dirname(host_removal_reference)
    fasta_file_name = os.path.basename(host_removal_reference)
    minimap2_index_cmd = f"singularity exec -B {hostremoval_reference_directory}:/mnt {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate host_removal && minimap2 -d {host_removal_out_folder}/{reference_prefix}.min /mnt/{fasta_file_name}\""
    try:
        print(f"Building minimap2 index for host genome {host_removal_reference}...")
        subprocess.run(minimap2_index_cmd, shell=True, check=True)
        print(f"Minimap2 index created: {host_removal_out_folder}/{reference_prefix}.min")
    except subprocess.CalledProcessError as e:
        print(f"Error during minimap2 index creation: {e}")
        return
    
    for fasta_file in fasta_files:
        file_prefix = os.path.splitext(os.path.basename(fasta_file))[0]

        minimap2_align_cmd = f"singularity exec {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate host_removal && minimap2 -ax map-ont -t {num_threads} {host_removal_out_folder}/{reference_prefix}.min {fasta_file} -o {host_removal_out_folder}/{file_prefix}_minimap.sam\""
        try:
            print(f"Aligning {fasta_file} to host genome...")
            subprocess.run(minimap2_align_cmd, shell=True, check=True)
            print(f"Alignment completed: {host_removal_out_folder}/{file_prefix}_minimap.sam")
        except subprocess.CalledProcessError as e:
            print(f"Error during alignment: {e}")
            return

        extract_unmaped_cmd = f"singularity exec -B {hostremoval_reference_directory}:/mnt {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate host_removal && samtools view -bS -@ {num_threads} -T /mnt/{fasta_file_name} -f 4 {host_removal_out_folder}/{file_prefix}_minimap.sam > {host_removal_out_folder}/{file_prefix}_unmaped_minimap.bam\""
        try:
            print(f"Extracting unmaped reads from {host_removal_out_folder}/{file_prefix}_minimap.sam...")
            subprocess.run(extract_unmaped_cmd, shell=True, check=True)
            print(f"Unmaped reads extracted: {host_removal_out_folder}/{file_prefix}_unmaped_minimap.bam")
        except subprocess.CalledProcessError as e:
            print(f"Error during unmaped reads extraction: {e}")
            continue

        sort_bam_cmd = f"singularity exec {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate host_removal && samtools sort -n {host_removal_out_folder}/{file_prefix}_unmaped_minimap.bam -o {host_removal_out_folder}/{file_prefix}_unmaped_sorted_minimap.bam\""
        try:
            print(f"Sorting BAM files...")
            subprocess.run(sort_bam_cmd, shell=True, check=True)
            print(f"Sorted BAM files: {host_removal_out_folder}/{file_prefix}_unmaped_sorted_minimap.bam")
        except subprocess.CalledProcessError as e:
            print(f"Error during BAM file sorting: {e}")
            continue
        
        bam_to_fastq_cmd = f"singularity exec {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate host_removal && bedtools bamtofastq -i {host_removal_out_folder}/{file_prefix}_unmaped_sorted_minimap.bam -fq {host_removal_out_folder}/{file_prefix}_fitted_raw.fastq\""
        try:
            print(f"Converting BAM files to FASTQ files...")
            subprocess.run(bam_to_fastq_cmd, shell=True, check=True)
            print(f"Conversion completed: {host_removal_out_folder}/{file_prefix}_fitted_raw.fastq")
        except subprocess.CalledProcessError as e:
            print(f"Error during BAM to FASTQ conversion: {e}")
            continue

def centrifuge(centrifuge_db, num_threads, sif_path):
    centrifuge_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "centrifuge_out")
    if not os.path.exists(centrifuge_out_folder):
        os.makedirs(centrifuge_out_folder)
        print(f"Created folder: {centrifuge_out_folder}")
    else:
        print(f"Folder already exists: {centrifuge_out_folder}")

    host_removal_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "host_removal_out")
        
    for fastq_file in os.listdir(host_removal_out_folder):
        if fastq_file.endswith("_fitted_raw.fastq"):
            file_prefix = os.path.splitext(fastq_file)[0]
            fastq_path = os.path.join(host_removal_out_folder, fastq_file)
            centrifuge_report_path = os.path.join(centrifuge_out_folder, f"{file_prefix.replace('_fitted_raw', '')}_report")
            centrifuge_classification_path = os.path.join(centrifuge_out_folder, f"{file_prefix.replace('_fitted_raw', '')}_result")
                   
            centrifuge_db_directory = os.path.dirname(centrifuge_db)
            db_name = os.path.basename(centrifuge_db)
            centrifuge_cmd = (
                f"singularity exec -B {centrifuge_db_directory}:/mnt {sif_path} /tools/Software/centrifuge-1.0.4/bin/centrifuge "
                f"-p {num_threads} -x /mnt/{db_name} "
                f"-q {fastq_path} "
                f"--report-file {centrifuge_report_path} "
                f"-S {centrifuge_classification_path}"
            )   
            try:
                print(f"Running Centrifuge for {fastq_file}...")
                subprocess.run(centrifuge_cmd, shell=True, check=True)
                print(f"Centrifuge classification completed for {fastq_file}")
            except subprocess.CalledProcessError as e:
                print(f"Error during Centrifuge classification for {fastq_file}: {e}")
                continue
    
            centrifuge_to_kraken_cmd = (
                f"singularity exec -B {centrifuge_db_directory}:/mnt {sif_path} /tools/Software/centrifuge-1.0.4/bin/centrifuge-kreport "
                f"-x /mnt/{db_name} "
                f"{centrifuge_classification_path} > "
                f"{centrifuge_out_folder}/{file_prefix}_kraken_report"
            )
            
            try:
                print(f"Converting Centrifuge result to Kraken format for {fastq_file}...")
                subprocess.run(centrifuge_to_kraken_cmd, shell=True, check=True)
                print(f"Centrifuge result converted to Kraken report for {fastq_file}")
            except subprocess.CalledProcessError as e:
                print(f"Error during Centrifuge-to-Kraken conversion for {fastq_file}: {e}")
                continue

def kraken2(kraken2_db, num_threads, sif_path):
    kraken2_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "kraken2_out")
    if not os.path.exists(kraken2_out_folder):
        os.makedirs(kraken2_out_folder)
        print(f"Created folder: {kraken2_out_folder}")
    else:
        print(f"Folder already exists: {kraken2_out_folder}")
    
    host_removal_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "host_removal_out")
    
    for fastq_file in os.listdir(host_removal_out_folder):
        if fastq_file.endswith("_fitted_raw.fastq"):
            file_prefix = os.path.splitext(fastq_file)[0]
            fastq_path = os.path.join(host_removal_out_folder, fastq_file)
            kraken2_report_path = os.path.join(kraken2_out_folder, f"{file_prefix.replace('_fitted_raw', '')}_kraken2_report")
            kraken2_classification_path = os.path.join(kraken2_out_folder, f"{file_prefix.replace('_fitted_raw', '')}_kraken2_result")
                   
            kraken2_cmd = (
                f"singularity exec -B {kraken2_db}:/mnt {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate base && /tools/Software/kraken2-2.1.3/kraken2 "
                f"--db /mnt "
                f"--report {kraken2_report_path} "
                f"--output {kraken2_classification_path} "
                f"--threads {num_threads} {fastq_path} \""
            )   
            try:
                print(f"Running Kraken2 for {fastq_file}...")
                subprocess.run(kraken2_cmd, shell=True, check=True)
                print(f"Kraken2 classification completed for {fastq_file}")
            except subprocess.CalledProcessError as e:
                print(f"Error during Kraken2 classification for {fastq_file}: {e}")
                continue

def abricate(num_threads, sif_path):
    def create_folder(folder):
        if not os.path.exists(folder):
            try:
                os.makedirs(folder)
                print(f"Created folder: {folder}")
            except OSError as e:
                print(f"Error creating folder {folder}: {e}")
                sys.exit(1)
        else:
            print(f"Folder already exists: {folder}")

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

    abricate_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "abricate_out")
    folders = {
        "arg": os.path.join(abricate_out_folder, "abricate_arg_out"),
        "vf": os.path.join(abricate_out_folder, "abricate_vf_out"),
        "is": os.path.join(abricate_out_folder, "abricate_is_out"),
        "abundance": os.path.join(abricate_out_folder, "arg_abundance_out")
    }
    for folder in folders.values():
        create_folder(folder)

    host_removal_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "host_removal_out")

    fasta_files = []
    for fastq_file in os.listdir(host_removal_out_folder):
        if fastq_file.endswith("_fitted_raw.fastq"):
            file_prefix = os.path.splitext(fastq_file)[0]
            fastq_path = os.path.join(host_removal_out_folder, fastq_file)
            fasta_output_path = os.path.join(abricate_out_folder, f"{file_prefix.replace('_fitted_raw', '')}.fasta")
            fastq_to_fasta_cmd = (
                f"singularity exec {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"seqkit fq2fa {fastq_path} > {fasta_output_path}\""
            )
            try:
                subprocess.run(fastq_to_fasta_cmd, shell=True, check=True)
                print(f"Converted {fastq_path} to {fasta_output_path}")
                fasta_files.append(fasta_output_path)
            except subprocess.CalledProcessError as e:
                print(f"Error converting {fastq_path} to {fasta_output_path}: {e}")
                continue

    genome_sizes = {}
    for fasta_file in fasta_files:
        file_prefix = get_file_prefix(os.path.basename(fasta_file), ['.fasta'])
        seqkit_cmd = (
            f"singularity exec {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
            f"seqkit stats {fasta_file} --all --tabular > {folders['abundance']}/{file_prefix}_genome.csv\""
        )
        try:
            subprocess.run(seqkit_cmd, shell=True, check=True)
            excel_file = os.path.join(folders['abundance'], f"{file_prefix}_genome.csv")
            genome_size = get_genome_size_from_file(excel_file)
            if genome_size is not None:
                genome_sizes[file_prefix] = genome_size
                print(f"Genome size of {fasta_file}: {genome_size} bp")
            else:
                print(f"Could not find genome size for {fasta_file}")
        except subprocess.CalledProcessError as e:
            print(f"Error getting genome size for {fasta_file}: {e}")
            continue

    for fasta_file in fasta_files:
        file_prefix = get_file_prefix(os.path.basename(fasta_file), ['.fasta'])

        abricate_arg_cmd = (
            f"singularity exec {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate abricate && abricate "
            f"--db ncbi --threads {num_threads} {fasta_file} > {folders['arg']}/{file_prefix}_ncbi_result\""          
        )
        try:
            print(f"Annotating ARGs of {fasta_file} ...")
            subprocess.run(abricate_arg_cmd, shell=True, check=True)
            print(f"Annotating ARGs completed: {folders['arg']}/{file_prefix}_ncbi_result")
        except subprocess.CalledProcessError as e:
            print(f"Error during annotating ARGs: {e}")
            continue
        
        abricate_vf_cmd = (
            f"singularity exec {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate abricate && abricate "
            f"--db vfdb --threads {num_threads} {fasta_file} > {folders['vf']}/{file_prefix}_vf_result\""
        )
        try:
            print(f"Annotating VFs of {fasta_file} ...")
            subprocess.run(abricate_vf_cmd, shell=True, check=True)
            print(f"Annotating VFs completed: {folders['vf']}/{file_prefix}_vf_result")
        except subprocess.CalledProcessError as e:
            print(f"Error during annotating VFs: {e}")
            continue

        abricate_is_cmd = (
            f"singularity exec {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate abricate && abricate "
            f"--db ISfinder --threads {num_threads} {fasta_file} > {folders['is']}/{file_prefix}_is_result\""
        )
        try:
            print(f"Annotating IS of {fasta_file} ...")
            subprocess.run(abricate_is_cmd, shell=True, check=True)
            print(f"Annotating IS completed: {folders['is']}/{file_prefix}_is_result")
        except subprocess.CalledProcessError as e:
            print(f"Error during annotating IS: {e}")
            continue

    for ncbi_result_file in os.listdir(folders["arg"]):
        if ncbi_result_file.endswith("_ncbi_result"):
            file_prefix = get_file_prefix(ncbi_result_file, ['_ncbi_result'])
            ncbi_result_path = os.path.join(folders["arg"], ncbi_result_file)
            abundance_output_path = os.path.join(folders["abundance"], f"{file_prefix}_abundance_result")
            
            genome_size_bp = genome_sizes.get(file_prefix, "unknown")
            genome_size_gb = bp_to_gb(genome_size_bp)
            
            if genome_size_gb is None:
                print(f"Genome size for {file_prefix} is invalid, skipping abundance calculation.")
                continue

            genome_size_gb_value = f"{genome_size_gb:.3f}" 

            abundance_calculate_cmd = (
                f"singularity exec {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate base && python "
                f"/tools/Software/abundance_calculate.py --i {ncbi_result_path} --data_size {genome_size_gb_value} -p {file_prefix} -o {folders['abundance']}\""
            )
            try:
                print(f"Calculating abundance for {ncbi_result_file} with genome size {genome_size_gb_value} Gb ...")
                subprocess.run(abundance_calculate_cmd, shell=True, check=True)
                print(f"Abundance calculation completed: {abundance_output_path}")
            except subprocess.CalledProcessError as e:
                print(f"Error calculating abundance for {ncbi_result_file}: {e}")
                
def metaflye(num_threads, sif_path):
    metaflye_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "metaflye_out")
    if not os.path.exists(metaflye_out_folder):
        os.makedirs(metaflye_out_folder)
        print(f"Created folder: {metaflye_out_folder}")
    else:
        print(f"Folder already exists: {metaflye_out_folder}")
    
    host_removal_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "host_removal_out")
        
    for fastq_file in os.listdir(host_removal_out_folder):
        if fastq_file.endswith("_fitted_raw.fastq"):
            file_prefix = fastq_file.replace("_fitted_raw.fastq", "")
            fastq_path = os.path.join(host_removal_out_folder, fastq_file)
            unique_fastq_path = os.path.join(host_removal_out_folder, f"{file_prefix}_unique.fastq")

            seqkit_rmdup_cmd = (
                f"singularity exec {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"seqkit rmdup {fastq_path} -o {unique_fastq_path}\""
            )
            try:
                print(f"Removing duplicates from {fastq_file}...")
                subprocess.run(seqkit_rmdup_cmd, shell=True, check=True)
                print(f"Duplicates removed. Output file: {unique_fastq_path}")
            except subprocess.CalledProcessError as e:
                print(f"Error during seqkit processing for {fastq_file}: {e}")
                continue

    for fastq_file in os.listdir(host_removal_out_folder):
        if fastq_file.endswith("_unique.fastq"):
            file_prefix = fastq_file.replace("_unique.fastq", "")
            fastq_path = os.path.join(host_removal_out_folder, fastq_file)
            metaflye_out_path = os.path.join(metaflye_out_folder, f"{file_prefix}_flye_out")

            metaflye_cmd = (
                f"singularity exec {sif_path} /tools/Software/Flye-2.9.2/bin/flye "
                f"--nano-raw {fastq_path} "
                f"--threads {num_threads} "
                f"--meta -o {metaflye_out_path} "
            )
            try:
                print(f"Running metaFlye for {fastq_file}...")
                subprocess.run(metaflye_cmd, shell=True, check=True)
                print(f"metaFlye assembly completed for {fastq_file}")
            except subprocess.CalledProcessError as e:
                print(f"Error during metaFlye assembly for {fastq_file}: {e}")
                continue

def mv_metaflye(metaflye_out_folder):
    for result_folder in os.listdir(metaflye_out_folder):
        if result_folder.endswith("_flye_out"):
            file_prefix = result_folder.replace("_flye_out", "")
            result_folder_path = os.path.join(metaflye_out_folder, result_folder)
            
            assembly_fasta_path = os.path.join(result_folder_path, "assembly.fasta")
            
            if os.path.exists(assembly_fasta_path):
                new_fasta_name = f"{file_prefix}.fasta"
                new_fasta_path = os.path.join(result_folder_path, new_fasta_name)
                
                try:
                    shutil.move(assembly_fasta_path, new_fasta_path)
                    print(f"Renamed {assembly_fasta_path} to {new_fasta_path}")
                except Exception as e:
                    print(f"Error renaming {assembly_fasta_path}: {e}")
            else:
                print(f"assembly.fasta not found in {result_folder_path}")

def nextpolish(num_threads, sif_path):
    nextpolish_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "nextpolish_out")
    if not os.path.exists(nextpolish_out_folder):
        os.makedirs(nextpolish_out_folder)
        print(f"Created folder: {nextpolish_out_folder}")
    else:
        print(f"Folder already exists: {nextpolish_out_folder}")
    
    host_removal_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "host_removal_out")
    
    metaflye_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "metaflye_out")
    
    for fastq_file in os.listdir(host_removal_out_folder):
        if fastq_file.endswith("_fitted_raw.fastq"):
            file_prefix = fastq_file.replace("_fitted_raw.fastq", "")
            fastq_path = os.path.join(host_removal_out_folder, fastq_file)
            lgs_out_path = os.path.join(nextpolish_out_folder, f"{file_prefix}.fofn")

            try:
                with open(lgs_out_path, 'w') as lgs_file:
                    lgs_file.write(fastq_path + "\n")
                print(f"fofn file created successfully for {fastq_file}. Output file: {lgs_out_path}")
            except Exception as e:
                print(f"Error during fofn file creation for {fastq_file}: {e}")
                continue

            cfg_file_path = os.path.join(nextpolish_out_folder, f"{file_prefix}.cfg")
            metaflye_prefix_path = os.path.join(metaflye_out_folder, f"{file_prefix}_flye_out")
            genome_path = os.path.join(metaflye_prefix_path, f"{file_prefix}.fasta")
            workdir_path = os.path.join(nextpolish_out_folder, f"{file_prefix}_nextpolish_out")
            cfg_content = f"""\

[General]
job_type = local
job_prefix = nextPolish
task = best
rewrite = yes
rerun = 3
parallel_jobs = 6
multithread_jobs = 5
genome = {genome_path} #组装结果文件
genome_size = auto
workdir = {workdir_path}
polish_options = -p {num_threads}

[lgs_option]
lgs_fofn = {lgs_out_path}
lgs_options = -min_read_len 1k -max_depth 100
lgs_minimap2_options = -x map-ont
"""

            try:
                with open(cfg_file_path, 'w') as cfg_file:
                    cfg_file.write(cfg_content)
                print(f"cfg file created successfully. Output file: {cfg_file_path}")
            except Exception as e:
                print(f"Error during cfg file creation for {lgs_out_path}: {e}")
          
            nextpolish_cmd = (
                f"singularity exec {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"/tools/Software/NextPolish/nextPolish {cfg_file_path}\""
            )
            
            try:
                print(f"Performing calibration procedures for {file_prefix}.fasta...")
                subprocess.run(nextpolish_cmd, shell=True, check=True)
                print(f"Calibration procedures completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error during calibration procedures for {file_prefix}.fasta: {e}")
                continue

def mv_nextpolish(nextpolish_out_folder):
    for result_folder in os.listdir(nextpolish_out_folder):
        if result_folder.endswith("_nextpolish_out"):
            file_prefix = result_folder.replace("_nextpolish_out", "")
            result_folder_path = os.path.join(nextpolish_out_folder, result_folder)
            
            calibration_fasta_path = os.path.join(result_folder_path, "genome.nextpolish.fasta")
            
            if os.path.exists(calibration_fasta_path):
                new_fasta_name = f"{file_prefix}_nextpolish.fasta"
                new_fasta_path = os.path.join(result_folder_path, new_fasta_name)
                
                try:
                    shutil.move(calibration_fasta_path, new_fasta_path)
                    print(f"Renamed {calibration_fasta_path} to {new_fasta_path}")
                except Exception as e:
                    print(f"Error renaming {calibration_fasta_path}: {e}")
            else:
                print(f"genome.nextpolish.fasta not found in {result_folder_path}")

def semi_bin(num_threads, sif_path):
    semi_bin_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "semi_bin_out")
    bam_out_folder = os.path.join(semi_bin_out_folder, "bam_out")
    
    folders = {
        "bam": bam_out_folder
    }
    
    for folder in folders.values():
        if not os.path.exists(folder):
            try:
                os.makedirs(folder)
                print(f"Created folder: {folder}")
            except OSError as e:
                print(f"Error creating folder {folder}: {e}")
                sys.exit(1)
        else:
            print(f"Folder already exists: {folder}")

    nextpolish_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "nextpolish_out")
    
    prefixes = []
    
    for folder_name in os.listdir(nextpolish_out_folder):
        folder_path = os.path.join(nextpolish_out_folder, folder_name)
        
        if os.path.isdir(folder_path):
            if folder_name.endswith("_nextpolish_out"):
                prefix = folder_name.rsplit("_nextpolish_out", 1)[0]
                prefixes.append((folder_path, prefix))
    
    for folder_path, prefix in prefixes:
        fasta_files = []
        
        for file_name in os.listdir(folder_path):
            if file_name == f"{prefix}_nextpolish.fasta":
                fasta_files.append(os.path.join(folder_path, file_name))
        
        for fasta_file in fasta_files:
            minimap2_index_cmd = (
                f"singularity exec {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"minimap2 -d catalogue.mmi {fasta_file} \""
            )
            try:
                print(f"Creating index file for {fasta_file}...")
                subprocess.run(minimap2_index_cmd, shell=True, check=True)
                print(f"Index procedures completed successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error during index procedures for {fasta_file}: {e}")
                continue

    host_removal_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "host_removal_out")
    
    for fastq_file in os.listdir(host_removal_out_folder):
        if fastq_file.endswith("_fitted_raw.fastq"):
            file_prefix = fastq_file.replace("_fitted_raw.fastq", "")
            fastq_path = os.path.join(host_removal_out_folder, fastq_file)
            bam_out_path = os.path.join(bam_out_folder, f"{file_prefix}.bam")
            
            minimap2_cmd = (
                f"singularity exec {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"minimap2 -t {num_threads} -N 5 -ax map-ont catalogue.mmi {fastq_path} | samtools view -F 3584 -b --threads 8 > {bam_out_path} \""
            )
            try:
                print(f"BAM file generated for {fastq_path}...")
                subprocess.run(minimap2_cmd, shell=True, check=True)
                print(f"BAM file generated successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error during BAM file generated procedures for {fastq_path}: {e}")
                continue
            
            sorted_bam_out_path = os.path.join(bam_out_folder, f"{file_prefix}.sorted.bam")
            samtools_cmd = (
                f"singularity exec {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate base && "
                f"samtools sort -@ 10 {bam_out_path} > {sorted_bam_out_path} \""
            )
            try:
                print(f"BAM file sorted for {file_prefix}.bam...")
                subprocess.run(samtools_cmd, shell=True, check=True)
                print(f"BAM file sorted successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error during BAM file sorted procedures for {file_prefix}.bam: {e}")
                continue
    
    for folder_path, prefix in prefixes:
        fasta_files = []
        
        for file_name in os.listdir(folder_path):
            if file_name == f"{prefix}_nextpolish.fasta":
                fasta_files.append(os.path.join(folder_path, file_name))
    
        for fasta_file in fasta_files: 
            semi_bin_cmd = (
                f"singularity exec {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate SemiBin && "
                f"SemiBin single_easy_bin -i {fasta_file}  --sequencing-type long_read -b {sorted_bam_out_path} -o {semi_bin_out_folder}/{file_prefix}_bin_out --environment global  \""
            )
            try:
                print(f"Perform bin for {fasta_file}...")
                subprocess.run(semi_bin_cmd, shell=True, check=True)
                print(f"Perform bin successfully.")
            except subprocess.CalledProcessError as e:
                print(f"Error during bin procedures for {fasta_file}: {e}")
                continue

def checkm2(checkm2_db, num_threads, sif_path):
    checkm2_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "checkm2_out")
    if not os.path.exists(checkm2_out_folder):
        os.makedirs(checkm2_out_folder)
        print(f"Created folder: {checkm2_out_folder}")
    else:
        print(f"Folder already exists: {checkm2_out_folder}")

    semi_bin_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "semi_bin_out")
   
    prefixes = []
    
    for folder_name in os.listdir(semi_bin_out_folder):
        folder_path = os.path.join(semi_bin_out_folder, folder_name)
        
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
                        f"singularity exec -B {checkm2_db_directory}:/mnt {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate checkm2 && "
                        f"checkm2 predict --database_path /mnt/{db_file_name} --threads {num_threads} --input {output_bins_folder}/* --output-directory {checkm2_out_folder}/{prefix}_checkm2_out\""
                    )
                    try:
                        print(f"Perform checkm2 for {folder_path}...")
                        subprocess.run(checkm2_cmd, shell=True, check=True)
                        print(f"Perform checkm2 successfully.")
                    except subprocess.CalledProcessError as e:
                        print(f"Error during checkm2 procedures for {folder_path}: {e}")
                        continue

def gtdbtk(gtdbtk_db, num_threads, sif_path):
    gtdbtk_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "gtdbtk_out")
    if not os.path.exists(gtdbtk_out_folder):
        os.makedirs(gtdbtk_out_folder)
        print(f"Created folder: {gtdbtk_out_folder}")
    else:
        print(f"Folder already exists: {gtdbtk_out_folder}")

    semi_bin_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "semi_bin_out")

    prefixes = []
    
    for folder_name in os.listdir(semi_bin_out_folder):
        folder_path = os.path.join(semi_bin_out_folder, folder_name)
        
        if os.path.isdir(folder_path) and folder_name.endswith("_bin_out"):
            prefix = folder_name.rsplit("_bin_out", 1)[0]
            prefixes.append((folder_path, prefix))

            for subfolder_name in os.listdir(folder_path):
                output_bins_folder = os.path.join(folder_path, subfolder_name)
                
                if os.path.isdir(output_bins_folder) and subfolder_name == "output_bins":
                    print(f"Processing {output_bins_folder} for prefix {prefix}")

                    gtdbtk_cmd = (
                        f"singularity exec -B {gtdbtk_db}:/mnt {sif_path} /bin/bash -c \"source /tools/Miniconda3/bin/activate gtdbtk-2.2.6 && export GTDBTK_DATA_PATH=/mnt && "
                        f"gtdbtk classify_wf --genome_dir {output_bins_folder} --extension fa --skip_ani_screen --out_dir {gtdbtk_out_folder}/{prefix}_gtdbtk_out\""
                    )
                    try:
                        print(f"Perform gtdbtk for {folder_path}...")
                        subprocess.run(gtdbtk_cmd, shell=True, check=True)
                        print(f"Perform gtdbtk successfully.")
                    except subprocess.CalledProcessError as e:
                        print(f"Error during gtdbtk procedures for {folder_path}: {e}")
                        continue
                   
def main():
    parser = argparse.ArgumentParser(
        description="Execute the EasyNanoMeta-pipeline for nanopore metagenomic data analysis."
    )
    parser.add_argument(
        "-f", "--folder",
        required=True,
        help="Absolute path of the folder to search for fastq/fq files."
    )
    parser.add_argument(
        "-sif", "--sif-path",
        type=str,
        default=os.path.join(os.getcwd(), "easynanometa.sif"),
        help="Absolute path to the easynanometa.sif file. If not provided, defaults to easynanometa.sif in the current directory."
    )
    parser.add_argument(
        "-t", "--threads",
        type=int,
        default=24,
        help="Number of threads to use for all operations (default: 24)."
    )
    parser.add_argument(
        "-host-removal-reference", "--host-removal-reference",
        required=True,
        help="Path to the reference host genome fasta file for host removal."
    )
    parser.add_argument(
        "-centrifuge-db", "--centrifuge-db",
        required=True,
        help="Path to centrifuge database path."
    )
    parser.add_argument(
        "-kraken2-db", "--kraken2-db",
        required=True,
        help="Path to Kraken2 database path."
    )
    parser.add_argument(
        "-checkm2-db", "--checkm2-database",
        type=str,
        help="Path to checkm2 database path."
    )
    parser.add_argument(
        "-gtdbtk-db", "--gtdbtk-database",
        type=str,
        help="Path to GTDB-Tk database path."
    )    
    
    args = parser.parse_args()
       
    if not any(arg.startswith('-') for arg in sys.argv[1:]):
        parser.print_help()
        sys.exit(1)

    output_folder, adapters_removal_out_folder = create_output_folder()

    fastq_files = find_fastq_files(args.folder)

    for file_prefix, file_path in fastq_files.items():
        adapters_removal(file_path, file_prefix, output_folder, args.threads, args.sif_path)
    
    host_removal(args.host_removal_reference, args.threads, args.sif_path)
    
    centrifuge(args.centrifuge_db, args.threads, args.sif_path)
    
    kraken2(args.kraken2_db, args.threads, args.sif_path)
    
    abricate(args.threads, args.sif_path)
    
    metaflye(args.threads, args.sif_path)
    
    metaflye_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "metaflye_out")
    mv_metaflye(metaflye_out_folder)   
    
    nextpolish(args.threads, args.sif_path)

    nextpolish_out_folder = os.path.join(os.getcwd(), "easynanometa_result", "nextpolish_out")    
    mv_nextpolish(nextpolish_out_folder)

    semi_bin(args.threads, args.sif_path)

    checkm2(args.checkm2_database, args.threads, args.sif_path)
    
    gtdbtk(args.gtdbtk_database, args.threads, args.sif_path)

if __name__ == "__main__":
    main()
