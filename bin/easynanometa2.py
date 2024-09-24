#!/usr/bin/env python

# easynanometa2.py

import argparse
from dep import run_flye, run_kraken2, run_abricate, run_adapters_removal, run_host_removal, run_centrifuge, run_arg_abundance, run_nextpolish, run_semi_bin, run_checkm2, run_gtdbtk  # Import functions from dep.py

# Main function to choose Flye or Kraken2 or abricate
def main():
    # Create a top-level parser showing flye, kraken2 and abricate options
    parser = argparse.ArgumentParser(prog='easynanometa2.py', description='Batch execute Flye, Kraken2 or abricate for multiple fastq files.')
    subparsers = parser.add_subparsers(dest='tool', help='Choose the tool to run.')

    # Create a subcommand parser for flye
    parser_flye = subparsers.add_parser('flye', help='Run Flye assembler')
    parser_flye.add_argument('-i', '--input_dir', required=True, help='Directory containing the input fastq files.')
    parser_flye.add_argument('-o', '--output_dir', required=True, help='Directory for storing the Flye output.')
    parser_flye.add_argument('-s', '--singularity_image', required=True, help='Path to the Singularity image (easynanometa1_2.sif).')
    parser_flye.add_argument('-f', '--flye_executable', default='/tools/Software/Flye-2.9.2/bin/flye', help='Path to Flye executable within the Singularity image.')
    parser_flye.add_argument('-t', '--threads', default=24, type=int, help='Number of threads to use.')

    # Create a subcommand parser for kraken2
    parser_kraken2 = subparsers.add_parser('kraken2', help='Run Kraken2 classifier')
    parser_kraken2.add_argument('-i', '--input_dir', required=True, help='Directory containing the input fastq files.')
    parser_kraken2.add_argument('-o', '--output_dir', required=True, help='Directory for storing the Kraken2 output.')
    parser_kraken2.add_argument('-s', '--singularity_image', required=True, help='Path to the Singularity image (easynanometa1_2.sif).')
    parser_kraken2.add_argument('-k', '--kraken2_executable', default='/tools/Software/kraken2-2.1.3/kraken2', help='Path to Kraken2 executable within the Singularity image.')
    parser_kraken2.add_argument('-db', '--kraken2_db', default='/backup/database/kraken2', help='Path to the Kraken2 database.')
    parser_kraken2.add_argument('-t', '--threads', default=24, type=int, help='Number of threads to use.')

    # Create a subcommand parser for abricate
    parser_abricate = subparsers.add_parser('abricate', help='Run abricate to identify functional genes.')
    parser_abricate.add_argument('-i', '--input_dir', required=True, help='Directory containing the input .fa files.')
    parser_abricate.add_argument('-o', '--output_dir', required=True, help='Directory for storing the Abricate output.')
    parser_abricate.add_argument('-s', '--singularity_image', required=True, help='Path to the Singularity image (easynanometa.sif).')
    parser_abricate.add_argument('-e', '--abricate_executable', default='abricate', help='Path to Abricate executable within the Singularity image.')
    parser_abricate.add_argument('-d', '--db', default='ncbi', help='The database to use for Abricate (default: ncbi).')

    # Create a subcommand parser for adapters-removal
    parser_adapters_removal = subparsers.add_parser('adapters-removal', help='Run porechop_abi to remove adapters.')
    parser_adapters_removal.add_argument('-i', '--input_dir', required=True, help='Directory containing the input fastq files.')
    parser_adapters_removal.add_argument('-o', '--output_dir', required=True, help='Directory for storing the adapters-removal output.')
    parser_adapters_removal.add_argument('-s', '--singularity_image', required=True, help='Path to the Singularity image (easynanometa.sif).')
    parser_adapters_removal.add_argument('-e', '--adapters_removal_executable', default='porechop_abi', help='Path to porechop_abi executable within the Singularity image.')
    parser_adapters_removal.add_argument('-t', '--threads', default=24, type=int, help='Number of threads to use.')

    # Create a subcommand parser for host-removal
    parser_host_removal = subparsers.add_parser('host-removal', help='Run host_removal to remove host genome.')
    parser_host_removal.add_argument('-i', '--input_dir', required=True, help='Directory containing the adapers removal output.')
    parser_host_removal.add_argument('-o', '--output_dir', required=True, help='Directory for storing the host removal output.')
    parser_host_removal.add_argument('-s', '--singularity_image', required=True, help='Path to the Singularity image (easynanometa.sif).')
    parser_host_removal.add_argument('-t', '--threads', default=24, type=int, help='Number of threads to use.')
    parser_host_removal.add_argument('-r', '--host-removal-reference', required=True, help='Path to the reference host genome fasta file for host removal.')

    # Create a subcommand parser for Centrifuge
    parser_centrifuge = subparsers.add_parser('centrifuge', help='Run Centrifuge classifier.')
    parser_centrifuge.add_argument('-i', '--input_dir', required=True, help='Directory containing the fillted fastq files.')
    parser_centrifuge.add_argument('-o', '--output_dir', required=True, help='Directory for storing the centrifuge output.')
    parser_centrifuge.add_argument('-s', '--singularity_image', required=True, help='Path to the Singularity image (easynanometa.sif).')
    parser_centrifuge.add_argument('-t', '--threads', default=24, type=int, help='Number of threads to use.')
    parser_centrifuge.add_argument('-db', '--centrifuge_db', required=True, help='Path to centrifuge database path.')
    parser_centrifuge.add_argument('-c', '--centrifuge_executable', default='/tools/Software/centrifuge-1.0.4/bin/centrifuge', help='Path to centrifuge executable within the Singularity image.')
    parser_centrifuge.add_argument('-k', '--kreport_executable', default='/tools/Software/centrifuge-1.0.4/bin/centrifuge-kreport', help='Path to centrifuge-kreport executable within the Singularity image.')


    # Create a subcommand parser for ARGs abundance calulated
    parser_arg_abundance = subparsers.add_parser('arg-abundance', help='Run calculating abundance for ARGs.')
    parser_arg_abundance.add_argument('-i', '--input_dir', required=True, help='Directory containing the abricate results of ncbi.')
    parser_arg_abundance.add_argument('-o', '--output_dir', required=True, help='Directory for storing the ARGs abundance caculated output.')
    parser_arg_abundance.add_argument('-s', '--singularity_image', required=True, help='Path to the Singularity image (easynanometa.sif).')
    parser_arg_abundance.add_argument('-fa', '--fasta_dir', required=True, help='Path to the directory containing fasta flies.')

    # Create a subcommand parser for NextPolish
    parser_nextpolish = subparsers.add_parser('nextpolish', help='Run calibration for flye results.')
    parser_nextpolish.add_argument('-i', '--input_dir', required=True, help='Directory containing host removal output.')
    parser_nextpolish.add_argument('-o', '--output_dir', required=True, help='Directory for storing the nextpolish calibration output.')
    parser_nextpolish.add_argument('-s', '--singularity_image', required=True, help='Path to the Singularity image (easynanometa.sif).')
    parser_nextpolish.add_argument('-t', '--threads', default=24, type=int, help='Number of threads to use.')
    parser_nextpolish.add_argument('-f', '--flye_dir', required=True, help='Path to the directory containing fasta flies.')

    # Create a subcommand parser for SemiBin
    parser_semi_bin = subparsers.add_parser('semibin', help='Run SemiBin for NextPolish results.')
    parser_semi_bin.add_argument('-i', '--input_dir', required=True, help='Directory containing NextPolish output.')
    parser_semi_bin.add_argument('-o', '--output_dir', required=True, help='Directory for storing the SemiBin output.')
    parser_semi_bin.add_argument('-s', '--singularity_image', required=True, help='Path to the Singularity image (easynanometa.sif).')
    parser_semi_bin.add_argument('-t', '--threads', default=24, type=int, help='Number of threads to use.')
    parser_semi_bin.add_argument('-host-removal', '--host_removal_dir', required=True, help='Path to the directory containing host removal output.')

    # Create a subcommand parser for Checkm2
    parser_checkm2 = subparsers.add_parser('checkm2', help='Run Checkm2 for SemiBin output.')
    parser_checkm2.add_argument('-i', '--input_dir', required=True, help='Directory containing SemiBin output.')
    parser_checkm2.add_argument('-o', '--output_dir', required=True, help='Directory for storing the Checkm2 output.')
    parser_checkm2.add_argument('-s', '--singularity_image', required=True, help='Path to the Singularity image (easynanometa.sif).')
    parser_checkm2.add_argument('-t', '--threads', default=24, type=int, help='Number of threads to use.')
    parser_checkm2.add_argument('-db', '--checkm2_db', required=True, help='Path to Checkm2 database path.')

    # Create a subcommand parser for GTDBTK
    parser_gtdbtk = subparsers.add_parser('gtdbtk', help='Run GTDBTK for SemiBin output.')
    parser_gtdbtk.add_argument('-i', '--input_dir', required=True, help='Directory containing SemiBin output.')
    parser_gtdbtk.add_argument('-o', '--output_dir', required=True, help='Directory for storing the GTDBTK output.')
    parser_gtdbtk.add_argument('-s', '--singularity_image', required=True, help='Path to the Singularity image (easynanometa.sif).')
    parser_gtdbtk.add_argument('-t', '--threads', default=24, type=int, help='Number of threads to use.')
    parser_gtdbtk.add_argument('-db', '--gtdbtk_db', required=True, help='Path to GTDBTK database path.')

    args = parser.parse_args()

    # Choose tool
    if args.tool == 'flye':
        run_flye(args.input_dir, args.output_dir, args.singularity_image, args.flye_executable, args.threads)
    elif args.tool == 'kraken2':
        run_kraken2(args.input_dir, args.output_dir, args.singularity_image, args.kraken2_executable, args.kraken2_db, args.threads)
    elif args.tool == 'abricate':
        run_abricate(args.input_dir, args.output_dir, args.singularity_image, args.abricate_executable, args.db)
    elif args.tool == 'adapters-removal':
        run_adapters_removal(args.input_dir, args.output_dir, args.singularity_image, args.adapters_removal_executable, args.threads)
    elif args.tool == 'host-removal':
        run_host_removal(args.host_removal_reference, args.input_dir, args.output_dir, args.singularity_image, args.threads)
    elif args.tool == 'centrifuge':
        run_centrifuge(args.centrifuge_executable, args.kreport_executable, args.centrifuge_db, args.input_dir, args.output_dir, args.singularity_image, args.threads)
    elif args.tool == 'arg-abundance':
        run_arg_abundance(args.input_dir, args.output_dir, args.singularity_image, args.fasta_dir)
    elif args.tool == 'nextpolish':
        run_nextpolish(args.input_dir, args.output_dir, args.singularity_image, args.threads, args.flye_dir)
    elif args.tool == 'semibin':
        run_semi_bin(args.input_dir, args.output_dir, args.singularity_image, args.threads, args.host_removal_dir)
    elif args.tool == 'checkm2':
        run_checkm2(args.input_dir, args.output_dir, args.singularity_image, args.threads, args.checkm2_db)
    elif args.tool == 'gtdbtk':
        run_gtdbtk(args.input_dir, args.output_dir, args.singularity_image, args.threads, args.gtdbtk_db)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
