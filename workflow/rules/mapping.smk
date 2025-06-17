import json

###############
# HPViewer
###############

rule hpv_viewer:
    priority: 10
    input:
        get_trimmed_fastq
    output:
        "hpv_viewer_repeatmasker/{sample}/temp/{sample}.bam",
        # temp("hpv_viewer_repeatmasker/{sample}/temp/{sample}.sam"),
        "hpv_viewer_repeatmasker/{sample}/{sample}_HPV_profile.txt",
        temp("hpv_viewer_repeatmasker/{sample}")
    params:
        samplename = "{sample}"
    shell:
        "python {config[HPViewer]} -m homology-mask -1 {input[0]} -2 {input[1]} -p {threads} -c 90 -o hpv_viewer_repeatmasker/{params.samplename}"


###############
# Checkpoint
###############


SAMPLES = (
    pd.read_csv(config["samples"], sep="\t")
    .set_index("sample_id", drop=False)
    .sort_index()
)

rule create_hpviewer_summary:
    input:
        metrics_files=expand("hpv_viewer_repeatmasker/{sample}/{sample}_HPV_profile.txt", sample=SAMPLES["sample_id"])
    output:
        outfile="hpv_viewer_repeatmasker/dom_genotype_summary.tsv"
    shell:
        """
        base_dir="hpv_viewer_repeatmasker"
        newfile="{output.outfile}"

        # Temporary file for initial aggregation
        temp_file=$(mktemp)

        # Loop through each folder in the base directory
        for folder in {input.metrics_files}; do
            # Extract the sample name from the folder path
            sample_name=$(basename $(dirname "$folder"))
            # Modify the sample name by removing "_HPV_profile" if present
            sample_name=$(echo "$sample_name" | sed 's/_HPV_profile//')

            # Check if the file exists
            if [ -f "$folder" ]; then
                # Check if the temp_file exists and has content
                if [ -s "$temp_file" ]; then
                    # File exists and has content, append without header
                    awk -v sample="$sample_name" 'BEGIN{{FS=OFS="\t"}} NR>1{{print sample, $0}}' "$folder" >> "$temp_file"
                else
                    # File does not exist or is empty, append with header
                    awk -v sample="$sample_name" 'BEGIN{{FS=OFS="\t"}} NR==1{{print "Sample", $0}} NR>1{{print sample, $0}}' "$folder" > "$temp_file"
                fi
            else
                echo "File not found: $folder"
            fi
        done

        # Process the temp_file to get the final output
        awk 'BEGIN{{FS=OFS="\t"}} NR==1{{header=$0; next}} 
             {{
                 if ($4 > 300) {{
                     if ($1 in max && $4 > max[$1]) {{
                         max[$1] = $4
                         line[$1] = $0
                     }} else if (!($1 in max)) {{
                         max[$1] = $4
                         line[$1] = $0
                     }}
                 }}
             }} 
             END{{ 
                 print header > "{output.outfile}"
                 for (sample in line) {{
                     print line[sample] >> "{output.outfile}"
                 }}
             }}' "$temp_file"

        # Clean up temporary file
        rm "$temp_file"
        """




checkpoint extract_hpviewer_summary:
    input:
        "hpv_viewer_repeatmasker/dom_genotype_summary.tsv"
    output:
        "hpv_viewer_repeatmasker/dom_genotype_summary.json"
    run:
        import json
        # Read the CSV file
        df = pd.read_csv(input[0],sep='\t')
        # Create a dictionary from column 1 and column 4
        samples_dict = dict(zip(df.iloc[:, 0], df.iloc[:, 1]))
        # Save dictionary to a JSON file
        with open(output[0], 'w') as file:
            json.dump(samples_dict, file)





def get_genotype_dir_from_sample(sample_id):
    # Path to the JSON file produced by the checkpoint
    json_file_path = checkpoints.extract_hpviewer_summary.get().output[0]
    
    # Read the JSON file
    with open(json_file_path, 'r') as file:
        samples_dict = json.load(file)
    
    # Extract the genotype for the given sample ID
    genotype = samples_dict.get(sample_id)
    dir = REF.loc[genotype][1]
    return dir

def get_genotype_dir_from_sample_hg(sample_id):
    # Path to the JSON file produced by the checkpoint
    json_file_path = checkpoints.extract_hpviewer_summary.get().output[0]
    
    # Read the JSON file
    with open(json_file_path, 'r') as file:
        samples_dict = json.load(file)
    
    # Extract the genotype for the given sample ID, default to "HPV16" if not found
    genotype = samples_dict.get(sample_id, "HPV16")
    
    # Get the directory from the REF DataFrame
    dir = REF.loc[genotype][1]
    return dir



################
## BWA alignment
################
# this is to get the human output
rule bwa_map:
    priority: 10
    input:
        # Assuming get_bwa_index and get_trimmed_fastq are functions returning file paths
        bwa_index=get_bwa_index(),
        fastq=get_trimmed_fastq,
        genotype_dir=lambda wildcards: get_sample_hpv_genotype(wildcards.sample)
    output:
        "raw_bam/{sample}.bam" #temp("raw_bam_virus/{sample}.bam")
    threads: 12
    log:
        "logs/{sample}_bwa_map_vir.log"
    shell:
        "(bwa mem -M -t {threads} {input.genotype_dir} {input.fastq} | "
        "samtools view -b -f 4 - | samtools collate -O - | samtools bam2fq - | "
        "bwa mem -M -p -t {threads} {input.bwa_index} - | " 
        "samtools view -b -f 2 -F 2828 --threads {threads} - | " # Removes all secondary alignments
        "samtools view -Sb --threads {threads} - | samtools sort - > {output} && "
        "samtools index {output}) 2> {log}"

REF = pd.read_csv(config["ref_files"], sep="\t", header = None, index_col = 0)

rule virus_bwa_map:
    priority: 10
    input:
        # Assuming get_bwa_index and get_trimmed_fastq are functions returning file paths
        bwa_index=get_bwa_index(),
        fastq=get_trimmed_fastq,
        genotype_dir=lambda wildcards: get_sample_hpv_genotype(wildcards.sample)
    output:
        "raw_bam_virus/{sample}.bam" #temp("raw_bam_virus/{sample}.bam")
    params:
        # Assuming get_genotype_dir_from_sample is a function
    threads: 12
    log:
        "logs/{sample}_bwa_map.log"
    shell:
        "(bwa mem -M -t {threads} {input.bwa_index} {input.fastq} | "
        "samtools view -b -f 4 - | samtools collate -O - | samtools bam2fq - | "
        "bwa mem -M -p -t {threads} {input.genotype_dir} - | " # Assuming get_genotype_dir_from_sample returns an index or directory
        "samtools view -b -f 2 -F 2828 --threads {threads} - | " # Removes all secondary alignments
        "samtools view -Sb --threads {threads} - | samtools sort - > {output} && "
        "samtools index {output}) 2> {log}"

# NO LONGER WILL USE
# rule virus_bwa_map_shifted:
#     input:
#         "raw_bam_virus/{sample}.bam",
#     output:
#         # temp("raw_bam/{sample}.bam")
#         discarded = temp("raw_bam_shifted/{sample}_discarded.bam"),
#         read_names = temp("raw_bam_shifted/{sample}_read_names.txt"),
#         temp1 = temp("raw_bam_shifted/{sample}_temp1.bam"),
#         kept = temp("raw_bam_shifted/{sample}_kept.bam"),
#         shifted = temp("raw_bam_shifted/{sample}_shifted.bam")
#     params:
#         # now these live under params, not input
#         genotype_dir_shifted = lambda wc: get_sample_hpv_genotype_shifted(wc.sample),
#         genotype_region         = lambda wc: get_sample_hpv_info_shifted(wc.sample)
#     threads: 12
#     log:
#         "logs/{sample}_bwa_map_shifted.log"
#     shell:
#         """
#         # Get all read names that overlap the region
#         samtools view {input} "{params.genotype_region}" | awk '{{print $1}}' | sort | uniq > {output.read_names}

#         # Extract both reads in each pair
#         samtools view -N {output.read_names} -b -h {input} -o {output.temp1}

#         # Separate kept and discarded reads
#         samtools view -N {output.read_names} -b {input} -U {output.discarded} -o {output.kept}

#         # bam to fastq conversion then bwa mem with the shifted genome
#         samtools collate -O {output.discarded} | samtools bam2fq - |
#         bwa mem -M -p -t {threads} {params.genotype_dir_shifted} - |
#         samtools view -b -f 2 -F 2828 --threads {threads} - | 
#         samtools view -Sb --threads {threads} - > {output.shifted}

#         """
    

        

##########################################
## raw bams without any filtering
## fixmate, sort, index and stats bam file
rule samtools_sort_index_stats:
    input:
        "raw_bam/{sample}.bam"
    output:
        #temp(bam = "raw_bam/{sample}_sorted.bam"), doesn't work
        #bai = "raw_bam/{sample}_sorted.bam.bai",
        #stat= "raw_bam/{sample}_sorted.bam.stats.txt"
        "raw_bam/{sample}_sorted.bam",
        "raw_bam/{sample}_sorted.bam.stats.txt"
    threads: 12
    shell:
        ## --threads flag failed 
        # "(samtools sort -n -@ {threads} -o - {input} |"
        # "samtools fixmate - - | "
        # "samtools sort  -@ {threads} -o {output[0]} && "
        # "samtools index -@ {threads} {output[0]} && "
        # "samtools stats -@ {threads} {output[0]} > {output[1]})"
        "samtools sort -n -@ {threads} -o - {input} |"
        "samtools fixmate - - | "
        "samtools sort -@ {threads} -o {output[0]} && "
        "samtools index -@ {threads} {output[0]} && "
        "samtools stats -@ {threads} {output[0]} > {output[1]}"


rule samtools_sort_index_stats_virus:
    input:
        "raw_bam_virus/{sample}.bam"
    output:
        #temp(bam = "raw_bam/{sample}_sorted.bam"), doesn't work
        #bai = "raw_bam/{sample}_sorted.bam.bai",
        #stat= "raw_bam/{sample}_sorted.bam.stats.txt"
        "raw_bam_virus/kept_{sample}_sorted.bam",
        "raw_bam_virus/kept_{sample}_sorted.bam.stats.txt"
    threads: 12
    shell:
        ## --threads flag failed
        # "(samtools fixmate -@ {threads} -m {input} - | "
        # "samtools sort -@ {threads} -o {output[0]} {input} && "
        # "samtools index -@ {threads} {output[0]} && "
        # "samtools stats -@ {threads} {output[0]} > {output[1]}"
        "(samtools sort -n -@ {threads} -o - {input} |"
        "samtools fixmate - - | "
        "samtools sort  -@ {threads} -o {output[0]} && "
        "samtools index -@ {threads} {output[0]} && "
        "samtools stats -@ {threads} {output[0]} > {output[1]})"


# rule samtools_sort_index_stats_shifted:
#     input:
#         "raw_bam_shifted/{sample}_shifted.bam"
#     output:
#         #temp(bam = "raw_bam/{sample}_sorted.bam"), doesn't work
#         #bai = "raw_bam/{sample}_sorted.bam.bai",
#         #stat= "raw_bam/{sample}_sorted.bam.stats.txt"
#         "raw_bam_shifted/{sample}_sorted.bam",
#         "raw_bam_shifted/{sample}_sorted.bam.stats.txt"
#     threads: 12
#     shell:
#         ## --threads flag failed
#         # "(samtools sort -n -@ {threads} -o - {input} |"
#         # "samtools fixmate - - | "
#         # "samtools sort  -@ {threads} -o {output[0]} && "
#         # "samtools index -@ {threads} {output[0]} && "
#         # "samtools stats -@ {threads} {output[0]} > {output[1]})"
#         "(samtools sort -n -@ {threads} -o - {input} |"
#         "samtools fixmate - - | "
#         "samtools sort  -@ {threads} -o {output[0]} && "
#         "samtools index -@ {threads} {output[0]} && "
#         "samtools stats -@ {threads} {output[0]} > {output[1]})"





##########################
# Consensus Cruncher
###########################

rule create_ini_file:
    input:
        "raw_bam/{sample}_sorted.bam"
    params:
        samplename = "{sample}"
    output:
        ini_file = "cc_ini/{sample}.ini"
    shell:
        """
        cat <<EOF > {output.ini_file}
        [fastq2bam]
        fastq1 = skip
        fastq2 = skip
        output = skip
        name = {params.samplename}
        bwa = /cluster/tools/software/bwa/0.7.15/bwa
        ref = bwa_ref
        samtools = /cluster/tools/software/rocky9/samtools/1.20/bin/samtools
        bpattern = NNT
        [consensus]
        bam = {input}
        c_output = cc_data
        """

# rule create_ini_file_shifted:
#     input:
#         "raw_bam_shifted/{sample}_sorted.bam"
#     params:
#         samplename = "{sample}"
#     output:
#         ini_file = "cc_ini/{sample}_shifted.ini"
#     shell:
#         """
#         cat <<EOF > {output.ini_file}
#         [fastq2bam]
#         fastq1 = skip
#         fastq2 = skip
#         output = skip
#         name = {params.samplename}
#         bwa = /cluster/tools/software/bwa/0.7.15/bwa
#         ref = bwa_ref
#         samtools = /cluster/tools/software/rocky9/samtools/1.20/bin/samtools
#         bpattern = NNT
#         [consensus]
#         bam = {input}
#         c_output = cc_data_shifted
#         """

rule create_ini_file_virus:
    input:
        "raw_bam_virus/kept_{sample}_sorted.bam"
    params:
        samplename = "kept_{sample}"
    output:
        ini_file = "cc_ini/{sample}_virus.ini"
    shell:
        """
        cat <<EOF > {output.ini_file}
        [fastq2bam]
        fastq1 = skip
        fastq2 = skip
        output = skip
        name = {params.samplename}
        bwa = /cluster/tools/software/bwa/0.7.15/bwa
        ref = bwa_ref
        samtools = /cluster/tools/software/rocky9/samtools/1.20/bin/samtools
        bpattern = NNT
        [consensus]
        bam = {input}
        c_output = cc_data_virus
        """



rule run_consensus_cruncher:
    input:
        ini_file = "cc_ini/{sample}.ini"
    output:
        consensus_output = "cc_data/{sample}_sorted/dcs_sc/{sample}_sorted.all.unique.dcs.sorted.bam"  # Update this as per your actual output file(s) or directory
    shell:
        """
        python3 {config[cc_run_hg]} -c {input.ini_file} consensus
        """

rule run_consensus_cruncher_virus:
    input:
        ini_file = "cc_ini/{sample}_virus.ini"
    output:
        consensus_output = "cc_data_virus/kept_{sample}_sorted/dcs_sc/kept_{sample}_sorted.all.unique.dcs.sorted.bam"  # Update this as per your actual output file(s) or directory
    shell:
        """
        python3 {config[cc_run_hpv]} -c {input.ini_file} consensus -b False
        """

# rule run_consensus_cruncher_shifted:
#     input:
#         ini_file = "cc_ini/{sample}_shifted.ini"
#     output:
#         consensus_output = "cc_data_shifted/{sample}_sorted/dcs_sc/{sample}_sorted.all.unique.dcs.sorted.bam"  # Update this as per your actual output file(s) or directory
#     shell:
#         """
#         python3 {config[cc_run_hpv]} -c {input.ini_file} consensus -b False
#         """

rule get_dedup_bam_from_cc:
    input:
        file = "cc_data/{sample}_sorted/dcs_sc/{sample}_sorted.all.unique.dcs.sorted.bam"
    output:
        dedup_bam = "dedup_bam_umi_pe/{sample}_dedup.bam",
    shell:
        """ 
        cp {input.file} {output.dedup_bam}
        cp {input.file}.bai {output.dedup_bam}.bai
        """


rule get_dedup_bam_from_cc_virus:
    input:
        #file = "cc_data_shifted/{sample}_sorted/dcs_sc/{sample}_sorted.all.unique.dcs.sorted.bam",
        file = "cc_data_virus/kept_{sample}_sorted/dcs_sc/kept_{sample}_sorted.all.unique.dcs.sorted.bam"
    output:
        dedup_bam = "dedup_bam_umi_pe_shifted/{sample}_dedup.bam",
        #dedup_bam_unsorted = temp("dedup_bam_umi_pe_shifted/{sample}_dedup_unsorted.bam")
    shell:
        """ 
        cp {input.file} {output.dedup_bam}
        cp {input.file}.bai {output.dedup_bam}.bai
        """

## 20240420 
# # editing this out for the time being but will reimplement when it's necessary to do.
# rule get_dedup_bam_from_cc_virus:
#     input:
#         file = "cc_data_kept/kept_{sample}_sorted/dcs_sc/kept_{sample}_sorted.all.unique.dcs.sorted.bam"
#     output:
#         dedup_bam = "dedup_bam_umi_pe_shifted/{sample}_dedup.bam",
#         #dedup_bam_unsorted = temp("dedup_bam_umi_pe_shifted/{sample}_dedup_unsorted.bam")
#     shell:
#         """ 
#         cp {input.file} {output.dedup_bam}
#         cp {input.file}.bai {output.dedup_bam}.bai
#         """

############################################
## infer insert size for paired-end reads_qc
############################################
rule insert_size:
    input:
        #"dedup_bam_pe/{sample}_dedup.bam"
        get_dedup_bam
    output:
        txt = "fragment_size/{sample}_insert_size_metrics.txt",
        hist = "fragment_size/{sample}_insert_size_histogram.pdf"
    params:
        pipeline_env = env_dir
    log:
        "logs/{sample}_picard_insert_size.log"
    shell:
         r"""
        set -euo pipefail
        module load R

        java -jar /cluster/tools/software/picard/2.10.9/picard.jar \
            CollectInsertSizeMetrics \
            M=0.01 \
            I={input} \
            O={output.txt} \
            H={output.hist} \
        2> {log}
        """

rule insert_size_virus:
    input:
        #"dedup_bam_pe/{sample}_dedup.bam"
        get_dedup_bam_virus
    output:
        txt = "fragment_size_virus/{sample}_insert_size_metrics.txt",
        hist = "fragment_size_virus/{sample}_insert_size_histogram.pdf"
    params:
        pipeline_env = env_dir
    log:
        "logs/{sample}_picard_insert_size_virus.log"
    shell:
         r"""
        set -euo pipefail
        module load R

        java -jar /cluster/tools/software/picard/2.10.9/picard.jar \
            CollectInsertSizeMetrics \
            M=0.01 \
            I={input} \
            O={output.txt} \
            H={output.hist} \
        2> {log}
        """



rule create_metrics_summary:
    input:
        metrics_files=
        expand("fragment_size/{samples}_insert_size_metrics.txt",samples = SAMPLES["sample_id"])
    output:
        outfile= "fragment_size/fragment_length_summary.csv"
    params:
        jobname = "custom_job_name_for_my_rule"
    shell:
        """
        # Create the header
        outfile={output.outfile}
        echo -n "filename,READ_PAIRS," > $outfile
        for i in $(seq 0 600); do
            echo -n "$i," >> $outfile
        done
        echo "" >> $outfile

        # Process each metrics file
        for file in {input.metrics_files}; do
            # Create an array with 601 zeroes
            arr=($(for i in $(seq 0 600); do echo -n "0 "; done))

            # Extract the filename
            filename=$(basename $file)

            # Extract the READ_PAIRS value
            read_pairs=$(awk '/READ_PAIRS/{{getline; print $7}}' $file)

            # Update the array with the actual data
            awk -v arr="${{arr[*]}}" -v filename="$filename" -v read_pairs="$read_pairs" \\
                'BEGIN{{split(arr, a, " ")}}
                $1 ~ /^[0-9]+$/ && $1 >= 0 && $1 <= 600 {{a[$1+1]=$2}}
                END{{printf("%s,%s,", filename, read_pairs); for(i=1; i<=601; i++) printf("%s,", a[i]); printf("\\n")}}' $file >> $outfile
        done
        """



rule create_metrics_summary_virus:
    input:
        metrics_files=expand("fragment_size_virus/{sample}_insert_size_metrics.txt", sample=SAMPLES["sample_id"])
    output:
        outfile= "fragment_size_virus/fragment_length_summary.csv"
    params:
        jobname = "custom_job_name_for_my_rule"
    shell:
        """
        # Create the header
        outfile={output.outfile}
        echo -n "filename,READ_PAIRS," > $outfile
        for i in $(seq 0 600); do
            echo -n "$i," >> $outfile
        done
        echo "" >> $outfile

        # Process each metrics file
        for file in {input.metrics_files}; do
            # Create an array with 601 zeroes
            arr=($(for i in $(seq 0 600); do echo -n "0 "; done))

            # Extract the filename
            filename=$(basename $file)

            # Extract the READ_PAIRS value
            read_pairs=$(awk '/READ_PAIRS/{{getline; print $7}}' $file)

            # Update the array with the actual data
            awk -v arr="${{arr[*]}}" -v filename="$filename" -v read_pairs="$read_pairs" \\
                'BEGIN{{split(arr, a, " ")}}
                $1 ~ /^[0-9]+$/ && $1 >= 0 && $1 <= 600 {{a[$1+1]=$2}}
                END{{printf("%s,%s,", filename, read_pairs); for(i=1; i<=601; i++) printf("%s,", a[i]); printf("\\n")}}' $file >> $outfile
        done
        """



######################################
# COVERAGE
#######################################

rule coverage_virus:
    input:
        #"dedup_bam_pe/{sample}_dedup.bam"
        get_dedup_bam_virus
    output:
        "HPV_coverage/{sample}_depths.txt",
    log:
        "logs/{sample}_samtools_depth.log"
    shell:
        """
        samtools depth -d 100000 {input} > {output}  
        """



########################
# DRAFT
##########################

rule filter_fragments_by_tlen:
    input:
        bam=get_dedup_bam_virus
    output:
        below150="filtered_bams/{sample}_below150bp.bam",
        above150="filtered_bams/{sample}_above150bp.bam",
        below150_index="filtered_bams/{sample}_below150bp.bam.bai",
        above150_index="filtered_bams/{sample}_above150bp.bam.bai"
    log:
        "logs/{sample}_filter_fragments.log"
    shell:
        r"""
        set -euo pipefail
        samtools view -h -f 2 {input.bam} | \
        awk 'BEGIN {{OFS="\t"}} 
             /^@/ {{print > "header.tmp"; next}} 
             ($9 <= -100 && $9 >= -150) || ($9 >= 100 && $9 <= 150) {{print > "below.tmp"}} 
             ($9 <= -151 && $9 >= -220) || ($9 >= 151 && $9 <= 220) {{print > "above.tmp"}}'

        cat header.tmp below.tmp | samtools view -b -o {output.below150}
        cat header.tmp above.tmp | samtools view -b -o {output.above150}

        rm header.tmp below.tmp above.tmp

        samtools index {output.below150}
        samtools index {output.above150}
        """


rule coverage_by_fragment_length:
    input:
        below150="filtered_bams/{sample}_below150bp.bam",
        above150="filtered_bams/{sample}_above150bp.bam"
    output:
        below150_depth="HPV_coverage/{sample}_below150bp_depths.txt",
        above150_depth="HPV_coverage/{sample}_above150bp_depths.txt"
    log:
        below150_log="logs/{sample}_below150bp_depth.log",
        above150_log="logs/{sample}_above150bp_depth.log"
    shell:
        r"""
        samtools depth -d 100000 {input.below150} > {output.below150_depth} 2> {log.below150_log}
        samtools depth -d 100000 {input.above150} > {output.above150_depth} 2> {log.above150_log}
        """





##################################
# END MOTIF
#####################################


rule end_motif_prep:
    input:
        sorted_bam =  "dedup_bam_umi_pe_shifted/{sample}_dedup.bam"
    output:
        bedpe = "end_motif/{sample}.bedpe"
    params:
        mapq=30
    resources: cpus=1, mem_mb=1500, time_min=60
    shell:
        """
        samtools sort -n {input.sorted_bam} -o temp_sorted_{wildcards.sample}.bam
        samtools view -bf 0x2 -q {params.mapq} temp_sorted_{wildcards.sample}.bam | \
        bedtools bamtobed -i stdin -bedpe > {output.bedpe}
        rm temp_sorted_{wildcards.sample}.bam
        """

rule run_rscript_motif_format:
    input:
        bedpe = "end_motif/{sample}.bedpe"
    output:
        fasta_5= temp( "end_motif/" + "{sample}_fasta_5.bed"),
        fasta_3= temp( "end_motif/" + "{sample}_fasta_3.bed")
    params:
        script = config["pipe_dir"] + "/scripts/motif_format_bedpe_hpv.R"
    resources: cpus=1, mem_mb=1500, time_min=60
    shell:
        "Rscript {params.script} "
        "--id {wildcards.sample} "
        "--bedpe {input.bedpe} "
        "--outdir end_motif"


rule compute_endMotif:
    input:
        bed3 = "end_motif/{sample}_fasta_3.bed",
        bed5 = "end_motif/{sample}_fasta_5.bed"
    output:
        fasta_out3 = "end_motif/{sample}_fasta_3_annotated.bed",
        fasta_out5 = "end_motif/{sample}_fasta_5_annotated.bed"
    params:
        genome_ref = lambda wildcards: get_sample_hpv_genotype(wildcards.sample)
    resources: cpus=1, mem_mb=1500, time_min=60
    shell:
        "bedtools getfasta -bedOut -fi {params.genome_ref} -bed {input.bed3} > {output.fasta_out3};"
        "bedtools getfasta -bedOut -fi {params.genome_ref} -bed {input.bed5} > {output.fasta_out5}"

rule run_rscript_motif_get_contexts:
    input:
        fasta_5=lambda wildcards:  "end_motif/" + wildcards.sample + "_fasta_5_annotated.bed",
        fasta_3=lambda wildcards:  "end_motif/" + wildcards.sample + "_fasta_3_annotated.bed"
    output:
        final_output =  "end_motif/{sample}_motifs.txt",
        final_output2 =  "end_motif/{sample}_raw.txt",
        final_output3 =  "end_motif/{sample}_MDS.txt"
    params:
        script = config["pipe_dir"] + "/scripts/motif_get_contexts.R"
    resources: cpus=1, mem_mb=1500, time_min=60
    shell:
        "Rscript {params.script} "
        "--id {wildcards.sample} "
        "--fasta_5 {input.fasta_5} "
        "--fasta_3 {input.fasta_3} "
        "--outdir end_motif"


rule MDS_summary:
    input:
        expand( "end_motif/{samples}_MDS.txt", samples=SAMPLES["sample_id"])
    output:
         "end_motif/MDS_scores.txt" 
    resources: cpus=1, mem_mb=1000, time_min=5
    shell:
        """
        for file in {input}; do
            sample=$(basename $file _MDS.txt)
            value=$(sed -n '2p' $file)
            echo "$sample $value" >> {output}
        done
        """




#####################################
# Coverage
#####################################

