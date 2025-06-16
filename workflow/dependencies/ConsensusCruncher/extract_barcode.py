#!/usr/bin/env python3

import os
import sys
import re
import argparse
import configparser
from subprocess import Popen, PIPE, call


def fastq2bam(args):
    """
    Extract molecular barcodes from paired-end sequencing reads using a barcode list,
    pattern, or the two combined. Remove constant spacer bases and combine paired
    barcodes before adding to the header of each read in FASTQ files.

    Barcode-extracted FASTQ files are written to the 'fastq_tag' directory and are
    subsequenntly aligned with BWA. BAM files are written to a 'bamfile' directory
    under the specified project folder.

    BARCODE DESIGN:
    You can input either a barcode list or barcode pattern or both. If both are provided, barcodes will first be matched
    with the list and then the constant spacer bases will be removed before the barcode is added to the header.

    N = random / barcode base
    A | C | G | T = constant spacer bases
    e.g. ATNNGT means barcode is flanked by two spacers matching 'AT' in front and 'GT' behind.
    """
    # Create directory for barcode extracted FASTQ files and BAM files
    fastq_dir = '{}/fastq_tag'.format(args.output)
    bam_dir = '{}/bamfiles'.format(args.output)

    # Check if dir exists and there's permission to write
    if not os.path.exists(fastq_dir) and os.access(args.output, os.W_OK):
        os.makedirs(fastq_dir)
    if not os.path.exists(bam_dir) and os.access(args.output, os.W_OK):
        os.makedirs(bam_dir)

    # Set file variables
    filename = os.path.basename(args.fastq1).split(args.name, 1)[0]
    outfile = "{}/{}".format(fastq_dir, filename)

    ####################
    # Extract barcodes #
    ####################
    if args.blist is not None and args.bpattern is not None:
        extractb_cmd = "{}/ConsensusCruncher/extract_barcodes.py --read1 {} --read2 {} --outfile {} --bpattern {} "
        "--blist {}".format(code_dir, args.fastq1, args.fastq2, outfile, args.bpattern, args.blist)
    elif args.blist is None:
        extractb_cmd = "{}/ConsensusCruncher/extract_barcodes.py --read1 {} --read2 {} --outfile {} --bpattern {}".format(
            code_dir, args.fastq1, args.fastq2, outfile, args.bpattern)
    else:
        extractb_cmd = "{}/ConsensusCruncher/extract_barcodes.py --read1 {} --read2 {} --outfile {} --blist {}".format(
            code_dir, args.fastq1, args.fastq2, outfile, args.blist)

    print(extractb_cmd)
    os.system(extractb_cmd)


    # Create directories for bad barcodes and barcode distribution histograms
    if args.blist is not None:
        bad_barcode_dir = '{}/fastq_tag/bad_barcode'.format(args.output)
        barcode_dist_dir = '{}/fastq_tag/barcode_dist'.format(args.output)

        if not os.path.exists(bad_barcode_dir) and os.access(args.output, os.W_OK):
            os.makedirs(bad_barcode_dir)

        if not os.path.exists(barcode_dist_dir) and os.access(args.output, os.W_OK):
            os.makedirs(barcode_dist_dir)

        # Move files 
        os.rename('{}/{}_r1_bad_barcodes.txt'.format(fastq_dir, filename),
              '{}/{}_r1_bad_barcodes.txt'.format(bad_barcode_dir, filename))
        os.rename('{}/{}_r2_bad_barcodes.txt'.format(fastq_dir, filename),
              '{}/{}_r2_bad_barcodes.txt'.format(bad_barcode_dir, filename))
        os.rename('{}/{}_barcode_stats.png'.format(fastq_dir, filename),
              '{}/{}_barcode_stats.png'.format(barcode_dist_dir, filename))


if __name__ == '__main__':
    # Set up mode parser (turn off help message, to be added later)
    main_p = argparse.ArgumentParser(add_help=False)
    main_p.add_argument('-c', '--config', default=None,
                        help="Specify config file. Commandline option overrides config file (Use config template).")

    # Parse out config file (sub_args) and other command line args (remaining_args) to override config
    sub_args, remaining_args = main_p.parse_known_args()

    # Re-initialize parser with help message enabled
    main_p = argparse.ArgumentParser(parents=[main_p], add_help=True,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    sub = main_p.add_subparsers(help='sub-command help', dest='subparser_name')

    # Mode help messages
    mode_fastq2bam_help = "Extract molecular barcodes from paired-end sequencing reads using a barcode list and/or " \
                          "a barcode pattern."
    mode_consensus_help = "Almalgamate duplicate reads in BAM files into single-strand consensus sequences (SSCS) and" \
                          " duplex consensus sequences (DCS). Single reads with complementary duplex strands can also" \
                          " be corrected with 'Singleton Correction'."

    # Add subparsers
    sub_a = sub.add_parser('fastq2bam', help=mode_fastq2bam_help)
    sub_b = sub.add_parser('consensus', help=mode_consensus_help)

    # fastq2bam arg help messages
    fastq1_help = "FASTQ containing Read 1 of paired-end reads. [MANDATORY]"
    fastq2_help = "FASTQ containing Read 2 of paired-end reads. [MANDATORY]"
    output_help = "Output directory, where barcode extracted FASTQ and BAM files will be placed in " \
                  "subdirectories 'fastq_tag' and 'bamfiles' respectively (dir will be created if they " \
                  "do not exist). [MANDATORY]"
    filename_help = "Output filename. If none provided, default will extract output name by taking everything left of" \
                    " '_R'."
    bwa_help = "Path to executable bwa. [MANDATORY]"
    samtools_help = "Path to executable samtools. [MANDATORY]"
    ref_help = "Reference (BWA index). [MANDATORY]"
    genome_help = "Genome version (e.g. hg19 or hg38), default: hg19"
    bpattern_help = "Barcode pattern (N = random barcode bases, A|C|G|T = fixed spacer bases). [MANDATORY]"
    blist_help = "List of barcodes (Text file with unique barcodes on each line). [MANDATORY]"
    bdelim_help = "Delimiter before barcode in read name " \
                  "(e.g. '|' in 'HWI-D00331:196:C900FANXX:7:1110:14056:43945|TTTT')"
    bam_help = "Bam files" #Jinfeng
    coutput_help = "coutput" #Jinfeng
    bedfile = "bed file" #Jinfeng
    cleanup_help = "clean up" #Jinfeng

    # Update subparsers with config
    if sub_args.config is not None:
        defaults = {"fastq1": fastq1_help,
                    "fastq2": fastq2_help,
                    "output": output_help,
                    "name": "_R",
                    "bwa": bwa_help,
                    "ref": ref_help,
                    "samtools": samtools_help,
                    "bpattern": None,
                    "blist": None,
                    "bam": bam_help,
                    "c_output": coutput_help,
                    "scorrect": 'True',
                    "genome": 'hg19',
                    "bedfile": bedfile,
                    "cutoff": 0.7,
                    "bdelim": '|',
                    "cleanup": cleanup_help}

        config = configparser.ConfigParser()
        config.read(sub_args.config)

        if config.has_section("fastq2bam"):
            # Add config file args to fastq2bam mode
            defaults.update(dict(config.items("fastq2bam")))
            sub_a.set_defaults(**defaults)
        if config.has_section("consensus"):
            # Add config file args to consensus mode
            defaults.update(dict(config.items("consensus")))
            sub_b.set_defaults(**defaults)

    # Parse commandline arguments
    sub_a.add_argument('--fastq1', dest='fastq1', metavar="FASTQ1", type=str, help=fastq1_help)
    sub_a.add_argument('--fastq2', dest='fastq2', metavar="FASTQ2", type=str, help=fastq2_help)
    sub_a.add_argument('-o', '--output', dest='output', type=str, help=output_help)
    sub_a.add_argument('-n', '--name', metavar="FILENAME", type=str, help=filename_help)
    sub_a.add_argument('-b', '--bwa', metavar="BWA", help=bwa_help, type=str)
    sub_a.add_argument('-r', '--ref', metavar="REF", help=ref_help, type=str)
    sub_a.add_argument('-s', '--samtools', metavar="SAMTOOLS", help=samtools_help, type=str)
    sub_a.add_argument('-p', '--bpattern', metavar="PATTERN", type=str, help=bpattern_help)
    sub_a.add_argument('-l', '--blist', metavar="LIST", type=str, help=blist_help)
    sub_a.set_defaults(func=fastq2bam)

