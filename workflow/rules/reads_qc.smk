###############################################
## get a copy of raw fq and rename to smaple ID
## or combined multiple lanes data and rename
## to make sure consistant wildcard.sample
################################################
## single end
# rule merge_and_rename_fq_se:
#     input:
#         get_raw_fastq_se
#     output:
#         temp("renamed_fq/{sample}.fastq.gz"),
#     shell:
#         "cat {input} > {output}"

## paired-end
rule merge_and_rename_fq_pe:
    input:
        R1 = get_raw_fastq_pe_R1,
        R2 = get_raw_fastq_pe_R2,
    output:
        # "renamed_fq/{sample}_R1.fastq.gz",
        # "renamed_fq/{sample}_R2.fastq.gz",
        temp("renamed_fq/{sample}_R1.fastq.gz"),
        temp("renamed_fq/{sample}_R2.fastq.gz"),
    shell:
        "cat {input.R1} > {output[0]} && "
        "cat {input.R2} > {output[1]} "

## Obsoleted !
umi_regex = r"""(?P<umi_1>\\^[ACGT]{6})"""
###########################################################
### extract UMI barcode and add it to FASTQ headers, p.r.n.
### pair-end unzipped FASTQ only with ConsensusCruncher
###########################################################
rule extract_barcode:
    input:
        get_renamed_fastq
    output:
        file1 = temp("barcoded_fq_pe/{sample}_unzip_R1.fastq"),
        file2 = temp("barcoded_fq_pe/{sample}_unzip_R2.fastq"),
        file3 = temp("barcoded_fq_pe/{sample}_R1.fastq.gz"),
        file4 = temp("barcoded_fq_pe/{sample}_R2.fastq.gz"),
        # file1 = "barcoded_fq_pe/{sample}_unzip_R1.fastq",
        # file2 = "barcoded_fq_pe/{sample}_unzip_R2.fastq",
        # file3 = "barcoded_fq_pe/{sample}_R1.fastq.gz",
        # file4 = "barcoded_fq_pe/{sample}_R2.fastq.gz",
    params:
        outfile = "barcoded_fq_pe/{sample}"          
    shell:
        """
        ## unzip gz files
        gunzip {input[0]} -c > {output.file1} && 
        gunzip {input[1]} -c > {output.file2}

        ## extract barcodes
        python3 {config[cc_extract]} \
        --blist /cluster/projects/scottgroup/people/jinfeng/data/barcode_DIME_5or6UMI.txt \
        --read1 {output.file1} --read2 {output.file2} \
        --outfile {params.outfile}

        ## gzip
        gzip {params.outfile}_barcode_R1.fastq && gzip {params.outfile}_barcode_R2.fastq
        # Optional: Rename the gzipped files if needed
        mv {params.outfile}_barcode_R1.fastq.gz {output.file3}
        mv {params.outfile}_barcode_R2.fastq.gz {output.file4}
        """




###########################################################
### extract UMI barcode and add it to FASTQ headers, p.r.n.
### UMI-tools can take care of single-end& pair-end!!
###########################################################
## paired-end
# rule umi_tools_extract_pe:
#     input:
#         get_renamed_fastq
#     output:
#         temp("barcoded_fq_pe/{sample}_R1.fastq.gz"),
#         temp("barcoded_fq_pe/{sample}_R2.fastq.gz"),
#         "barcoded_fq_pe/{sample}_extract.log"
#     params:
#         bcp = lambda wildcards: config["umi_pattern"]        ##  deactivate automatic wildcard expansion of {}
#     shell:
#         "umi_tools extract --extract-method=regex --stdin={input[0]} --read2-in={input[1]} "
#         "--bc-pattern={params.bcp} --bc-pattern2={params.bcp} "
#         "--stdout={output[0]} --read2-out={output[1]} --log={output[2]}"

# rule cc_extract_pe:
#     input:
#         get_renamed_fastq
#     output:
#         barcode_r1_fastq = "barcoded_fq_pe/{sample}/fastq_tag/{sample}_R1.fastq.gz",
#         barcode_r2_fastq = "barcoded_fq_pe/{sample}/fastq_tag/{sample}_R2.fastq.gz",
#         updated_stats = "barcoded_fq_pe/{sample}/fastq_tag/{sample}_extract.log"
#     params:
#         blist = get_blacklist(),
#         outfile = "barcoded_fq_pe/{sample}"
#     shell:
#         """
#         python3 /cluster/home/t116306uhn/Reference/ConsensusCruncher/ConsensusCruncher/extract_barcodes.py --read1 {input[0]} --read2 {input[1]} --outfile {params.outfile} --blist {params.blist}
        
#         qnameKey=$(head {input[0]} -n 1 | awk -F':' '{{print $1}}')
#         grep "^$qnameKey" {output.barcode_r1_fastq} | wc -l | awk '{{print "filter T R1: "$0}}' >> {output.updated_stats}
#         grep "^$qnameKey" {output.barcode_r2_fastq} | wc -l | awk '{{print "filter T R2: "$0}}' >> {output.updated_stats}
#         """


###################################
### automatically trimming adapters
### -q 20  Quality gz_trimming
### mimimal length : 20 nt
###################################


#for paired-end reads
rule trim_galore_pe:
    input:
        get_fastq_4trim
    output:
        # "trimmed_fq/{sample}_R1_val_1.fq.gz", #SHOULD BE temp
        # "trimmed_fq/{sample}_R2_val_2.fq.gz", #SHOULD BE temp
        temp("trimmed_fq/{sample}_R1_val_1.fq.gz"), #SHOULD BE temp
        temp("trimmed_fq/{sample}_R2_val_2.fq.gz"), #SHOULD BE temp
        "trimmed_fq/{sample}_R1.fastq.gz_trimming_report.txt",
        "trimmed_fq/{sample}_R2.fastq.gz_trimming_report.txt",
    params:
        ## path needs to be full path: failed to recognize space
        path = work_dir + "/trimmed_fq"
    threads: 12
    log:
        "logs/{sample}_trim_galore_pe.log"
    shell:
        "(trim_galore -q 20 --stringency 3 --length 20 "
        "--cores {threads} --paired -o {params.path} {input}) 2> {log}"


##########################################
### FASTQC for raw and trimmed fastq reads
##########################################
#for single-end reads
rule fastqc_se:
    input:
        get_renamed_fastq,
        get_trimmed_fastq
    output:
        "fastqc_se/{sample}_fastqc.zip",
        "fastqc_se/{sample}_trimmed_fastqc.zip"
    run:
        for fq in input:
            shell("fastqc {} -t 8 --outdir fastqc_se/".format(fq))

#for paired-end reads
rule fastqc_pe:
    input:
        get_renamed_fastq,
        get_trimmed_fastq
    output:
        "fastqc_pe/{sample}_R1_fastqc.zip",
        "fastqc_pe/{sample}_R2_fastqc.zip",
        "fastqc_pe/{sample}_R1_val_1_fastqc.zip",
        "fastqc_pe/{sample}_R2_val_2_fastqc.zip"
    run:
        for fq in input:
            shell("fastqc {} -t 8 --outdir fastqc_pe/".format(fq))
