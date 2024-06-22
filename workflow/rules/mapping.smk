import json

################
## BWA alignment
################
# this is to get the human output
rule bwa_map:
    priority: 10
    input:
        #config["bwa_index"],
        get_bwa_index(),
        get_trimmed_fastq,
        genotype_dir=lambda wildcards: get_genotype_dir_from_sample_hg(wildcards.sample)
    output:
        temp("raw_bam/{sample}.bam")
        # "raw_bam/{sample}.bam"
    threads: 12
    log:
        "logs/{sample}_bwa_map.log"
    shell:
        "(bwa mem -M -t {threads} {input.genotype_dir} {input[1]} {input[2]} | "
        "samtools view -b -f 4 - | samtools collate -O - | samtools bam2fq - | "
        "bwa mem -M -p -t {threads} {input[0]} - | "
        "samtools view -b -f 2 -F 2828 --threads {threads} - | " #this is removing all the secondary alignment, just keep in mind
        #" samtools view -b -f 4 - | samtools sort - | samtools bam2fq - | "
        # "bwa mem -M -p -t {threads} /cluster/projects/scottgroup/people/jinfeng/HPV-seq/bwa_HPVs/HPV16.fasta - |"
        "samtools view -Sb --threads {threads} - > {output}) 2> {log}"

REF = pd.read_csv(config["ref_files"], sep="\t", header = None, index_col = 0)

def get_sample_ids_from_checkpoint(wildcards):
    # Wait for the checkpoint to be ready and get its output
    hpvviewer_json = checkpoints.extract_hpviewer_summary.get(**wildcards).output[0]
    with open(hpvviewer_json, 'r') as file:
        samples_dict = json.load(file)
    # Extract the sample IDs
    sample_ids = list(samples_dict.keys())
    return sample_ids

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

rule virus_bwa_map:
    priority: 10
    input:
        # Assuming get_bwa_index and get_trimmed_fastq are functions returning file paths
        bwa_index=get_bwa_index(),
        fastq=get_trimmed_fastq,
        genotype_dir=lambda wildcards: get_genotype_dir_from_sample(wildcards.sample)
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

rule bwa_map_fusion:
    priority: 10
    input:
        #config["bwa_index"],
        get_bwa_fusion_index(),
        get_trimmed_fastq
    output:
        # temp("raw_bam/{sample}.bam")
        "raw_bam_fusion/{sample}.bam"
    threads: 12
    log:
        "logs/{sample}_bwa_map.log"
    shell:
        "(bwa mem -M -t {threads} {input} | "
        "samtools view -b -f 2 -F 524 --threads {threads} - | " #this is removing all the secondary alignment, just keep in mind
        #" samtools view -b -f 4 - | samtools sort - | samtools bam2fq - | "
        # "bwa mem -M -p -t {threads} /cluster/projects/scottgroup/people/jinfeng/HPV-seq/bwa_HPVs/HPV16.fasta - |"
        "samtools view -Sb --threads {threads} - > {output}) 2> {log}"

rule bwa_map_shifted:
    input:
        "raw_bam_virus/{sample}.bam"
    output:
        # temp("raw_bam/{sample}.bam")
        discarded = temp("raw_bam_shifted/{sample}_discarded.bam"),
        read_names = temp("raw_bam_shifted/{sample}_read_names.txt"),
        temp1 = temp("raw_bam_shifted/{sample}_temp1.bam"),
        kept = temp("raw_bam_shifted/{sample}_kept.bam"),
        shifted = temp("raw_bam_shifted/{sample}_shifted.bam")
    threads: 12
    log:
        "logs/{sample}_bwa_map_shifted.log"
    shell:
        """
        samtools view -b -h {input} "gi|333031|gb|K02718.1|PPH16:1000-6000" > {output.temp1}

        # Extract the read names from input.bam and write them to read_names.txt
        samtools view {output.temp1} | awk '{{print $1}}' > {output.read_names}

        # Filter the reads in input.bam using the read names in read_names.txt, write unselected reads to discarded.bam, and output the selected reads to filtered_output.bam
        samtools view -N {output.read_names} -b {input} -U {output.discarded} -o {output.kept}
        
        # bam to fastq conversion then bwa mem with the shifted genome
        samtools collate -O {output.discarded} | samtools bam2fq - |
        bwa mem -M -p -t {threads} /cluster/home/t116306uhn/Reference/HPV_shifted_genome_4000bp/shifted.fasta - |
        samtools view -b -f 2 -F 2828 --threads {threads} - | 
        samtools view -Sb --threads {threads} - > {output.shifted}
        """
    


rule hpv_viewer_repeatmasker:
    priority: 10
    input:
        get_trimmed_fastq
    output:
        "hpv_viewer_repeatmasker/{sample}/temp/{sample}.bam",
        # temp("hpv_viewer_repeatmasker/{sample}/temp/{sample}.sam"),
        "hpv_viewer_repeatmasker/{sample}/{sample}_HPV_profile.txt",
    params:
        samplename = "{sample}"
    shell:
        "python {config[HPViewer]} -m homology-mask -1 {input[0]} -2 {input[1]} -p {threads} -c 90 -o hpv_viewer_repeatmasker/{params.samplename}"


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
        "samtools fixmate {input} - | "
        "samtools sort -@ {threads} -o {output[0]} && "
        "samtools index -@ {threads} {output[0]} && "
        "samtools stats -@ {threads} {output[0]} > {output[1]}"

rule samtools_sort_index_stats_virus_kept:
    input:
        "raw_bam_shifted/{sample}_kept.bam"
    output:
        #temp(bam = "raw_bam/{sample}_sorted.bam"), doesn't work
        #bai = "raw_bam/{sample}_sorted.bam.bai",
        #stat= "raw_bam/{sample}_sorted.bam.stats.txt"
        "raw_bam_shifted/kept_{sample}_sorted.bam",
        "raw_bam_shifted/kept_{sample}_sorted.bam.stats.txt"
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


rule samtools_sort_index_stats_shifted:
    input:
        "raw_bam_shifted/{sample}_shifted.bam"
    output:
        #temp(bam = "raw_bam/{sample}_sorted.bam"), doesn't work
        #bai = "raw_bam/{sample}_sorted.bam.bai",
        #stat= "raw_bam/{sample}_sorted.bam.stats.txt"
        "raw_bam_shifted/{sample}_sorted.bam",
        "raw_bam_shifted/{sample}_sorted.bam.stats.txt"
    threads: 12
    shell:
        ## --threads flag failed
        # "(samtools sort -n -@ {threads} -o - {input} |"
        # "samtools fixmate - - | "
        # "samtools sort  -@ {threads} -o {output[0]} && "
        # "samtools index -@ {threads} {output[0]} && "
        # "samtools stats -@ {threads} {output[0]} > {output[1]})"
        "(samtools sort -n -@ {threads} -o - {input} |"
        "samtools fixmate - - | "
        "samtools sort  -@ {threads} -o {output[0]} && "
        "samtools index -@ {threads} {output[0]} && "
        "samtools stats -@ {threads} {output[0]} > {output[1]})"




##########################################################################
## to filter out unmapped & non-uniquely mapped, not properly paired reads
## Deduplication with markup, index and stats deduplicated file
###########################################################################

###############
## without UMIs
rule samtools_markdup_stats_pe:
    input:
        "raw_bam/{sample}_sorted.bam"
    output:
        bam = "dedup_bam_pe/{sample}_dedup.bam",
        #bai = "dedup_bam/{sample}_dedup.bam.bai",
        stat= "dedup_bam_pe/{sample}_dedup.bam.stats.txt"
    threads: 12
    shell:
        "(samtools view -b -f 2 -F 2828 --threads {threads} {input} | "
        "samtools markdup -@ {threads} -r - {output.bam} && "
        "samtools index -@ {threads} {output.bam} && "
        "samtools stats -@ {threads} {output.bam} > {output.stat})"

## single-end filtering differ
rule samtools_markdup_stats_se:
    input:
        "raw_bam/{sample}_sorted.bam"
    output:
        bam = "dedup_bam_se/{sample}_dedup.bam",
        #bai = "dedup_bam/{sample}_dedup.bam.bai",
        stat= "dedup_bam_se/{sample}_dedup.bam.stats.txt"
    threads: 12
    shell:
        "(samtools view -b -F 2820 --threads {threads} {input} | "
        "samtools markdup -@ {threads} -r - {output.bam} && "
        "samtools index -@ {threads} {output.bam} && "
        "samtools stats -@ {threads} {output.bam} > {output.stat})"

#######################################################################################
## with UMIs !!!!
## Deduplication with UMI-tools, which takes both UMI and coordinates info into account
## UMI-tools dosen't support parallel threads yet!!
## paired-end
# rule samtools_umi_tools_pe:
#     input:
#         "raw_bam/{sample}_sorted.bam"
#     output:
#         dedup_bam = "dedup_bam_umi_pe/{sample}_dedup.bam",
#         bam_stat = "dedup_bam_umi_pe/{sample}_dedup.bam.stats.txt",
#     params:
#         tmp_bam = "dedup_bam_umi_pe/{sample}_tmp.bam",
#         stat_prefix = "dedup_bam_umi_pe/{sample}_dedup"
#     threads: 12
#     log:
#         "logs/{sample}_dedup_umi.log"
#     shell:
#         ## --umi-separator='_' by default, could also be ":"
#         ## umi tools --method='unique', by default is 'directional'
#         ##  -output-stats can slow down the processing and increase memory usage considerably
#         #I modifed this line to filter out unmapped reads for only the HPV16 genome
#         "(umi_tools dedup --paired -I {input} -S {params.tmp_bam} --umi-separator='_' --output-stats={params.stat_prefix} && "
#         "samtools view -b -f 2 -F 2828 --threads {threads} {params.tmp_bam} > {output.dedup_bam} && "
#         # "samtools view -b -f 256 -F 2564 --threads {threads} {params.tmp_bam} > {output.dedup_bam} &&" #I added this
#         "samtools index -@ {threads} {output.dedup_bam}  && rm {params.tmp_bam} && "
#         "samtools stats -@ {threads} {output.dedup_bam} > {output.bam_stat}) 2> {log}"

rule samtools_secondary_umi_tools_pe:
    input:
        "raw_bam/{sample}_sorted.bam"
    output:
        dedup_bam = "dedup_secondary_bam_umi_pe/{sample}_dedup.bam",
        bam_stat = "dedup_secondary_bam_umi_pe/{sample}_dedup.bam.stats.txt",
    params:
        tmp_bam = "dedup_secondary_bam_umi_pe/{sample}_tmp.bam",
        stat_prefix = "dedup_secondary_bam_umi_pe/{sample}_dedup"
    threads: 12
    log:
        "logs/{sample}_dedup_umi.log"
    shell:
        ## --umi-separator='_' by default, could also be ":"
        ## umi tools --method='unique', by default is 'directional'
        ##  -output-stats can slow down the processing and increase memory usage considerably
        #I modifed this line to filter out unmapped reads for only the HPV16 genome
        "(umi_tools dedup --paired -I {input} -S {params.tmp_bam} --umi-separator='_' --output-stats={params.stat_prefix} && "
        "samtools view -b -f 257 -F 2572 --threads {threads} {params.tmp_bam} > {output.dedup_bam} &&" #I added this
        "samtools view -b -f 256 -F 2564 --threads {threads} {params.tmp_bam} > {output.dedup_bam} &&" #I added this
        "samtools index -@ {threads} {output.dedup_bam}  && rm {params.tmp_bam} && "
        "samtools stats -@ {threads} {output.dedup_bam} > {output.bam_stat}) 2> {log}"

## single-end with UMIs: different samtools filtering flags
rule samtools_umi_tools_se:
    input:
        "raw_bam/{sample}_sorted.bam"
    output:
        dedup_bam = "dedup_bam_umi_se/{sample}_dedup.bam",
        bam_stat = "dedup_bam_umi_se/{sample}_dedup.bam.stats.txt",
    params:
        tmp_bam = "dedup_bam_umi_se/{sample}_tmp.bam",
        stat_prefix = "dedup_bam_umi_se/{sample}_dedup"
    threads: 12
    log:
        "logs/{sample}_dedup_umi.log"
    shell:
        ## --umi-separator='_' by default, could also be ":"
        "(umi_tools dedup -I {input} -S {params.tmp_bam} --umi-separator=':' --output-stats={params.stat_prefix} && "
        "samtools view -b -F 2820 --threads {threads} {params.tmp_bam} > {output.dedup_bam} && "
        # "samtools view -b -f 256 -F 2564 --threads {threads} {params.tmp_bam} > {output.dedup_bam} &&" #I added this
        "samtools index -@ {threads} {output.dedup_bam}  && rm {params.tmp_bam} && "
        "samtools stats -@ {threads} {output.dedup_bam} > {output.bam_stat}) 2> {log}"


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
        samtools = /cluster/tools/software/samtools/1.9/bin/samtools
        bpattern = NNT
        [consensus]
        bam = {input}
        c_output = cc_data
        """

rule create_ini_file_shifted:
    input:
        "raw_bam_shifted/{sample}_sorted.bam"
    params:
        samplename = "{sample}"
    output:
        ini_file = "cc_ini/{sample}_shifted.ini"
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
        samtools = /cluster/tools/software/samtools/1.9/bin/samtools
        bpattern = NNT
        [consensus]
        bam = {input}
        c_output = cc_data_shifted
        """

rule create_ini_file_kept:
    input:
        "raw_bam_shifted/kept_{sample}_sorted.bam"
    params:
        samplename = "kept_{sample}"
    output:
        ini_file = "cc_ini/{sample}_kept_shifted.ini"
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
        samtools = /cluster/tools/software/samtools/1.9/bin/samtools
        bpattern = NNT
        [consensus]
        bam = {input}
        c_output = cc_data_kept
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

rule run_consensus_cruncher_kept:
    input:
        ini_file = "cc_ini/{sample}_kept_shifted.ini"
    output:
        consensus_output = "cc_data_kept/kept_{sample}_sorted/dcs_sc/kept_{sample}_sorted.all.unique.dcs.sorted.bam"  # Update this as per your actual output file(s) or directory
    shell:
        """
        python3 {config[cc_run_hpv]} -c {input.ini_file} consensus -b False
        """

rule run_consensus_cruncher_shifted:
    input:
        ini_file = "cc_ini/{sample}_shifted.ini"
    output:
        consensus_output = "cc_data_shifted/{sample}_sorted/dcs_sc/{sample}_sorted.all.unique.dcs.sorted.bam"  # Update this as per your actual output file(s) or directory
    shell:
        """
        python3 {config[cc_run_hpv]} -c {input.ini_file} consensus -b False
        """

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

## 20240618 
# editing this out for the time being but will reimplement when it's necessary to do.
# rule get_dedup_bam_from_cc_virus:
#     input:
#         file = "cc_data_shifted/{sample}_sorted/dcs_sc/{sample}_sorted.all.unique.dcs.sorted.bam",
#         file_kept = "cc_data_kept/kept_{sample}_sorted/dcs_sc/kept_{sample}_sorted.all.unique.dcs.sorted.bam"
#     output:
#         dedup_bam = "dedup_bam_umi_pe_shifted/{sample}_dedup.bam",
#         #dedup_bam_unsorted = temp("dedup_bam_umi_pe_shifted/{sample}_dedup_unsorted.bam")
#     shell:
#         """ 
#         samtools merge -f - {input.file} {input.file_kept} | \
#         samtools view -b -f 2 -F 2828 - | \
#         samtools sort -o {output.dedup_bam} -

#         samtools index {output.dedup_bam}
#         """

##
#this version just doesn't merge the samples from shifted and kept
rule get_dedup_bam_from_cc_virus:
    input:
        file = "cc_data_kept/kept_{sample}_sorted/dcs_sc/kept_{sample}_sorted.all.unique.dcs.sorted.bam"
    output:
        dedup_bam = "dedup_bam_umi_pe_shifted/{sample}_dedup.bam",
        #dedup_bam_unsorted = temp("dedup_bam_umi_pe_shifted/{sample}_dedup_unsorted.bam")
    shell:
        """ 
        cp {input.file} {output.dedup_bam}
        cp {input.file}.bai {output.dedup_bam}.bai
        """
############################################
## extract spike-ins bam after deduplication
############################################
## paired-end only so far !!
# rule samtools_spikein_sort_index_stats:
#     input:
#         #"raw_bam/{sample}_sorted.bam"      ## lead to ambiguous wildcards!?
#         #"dedup_bam_pe/{sample}_dedup.bam"
#         get_dedup_bam
#     output:
#         bam = "dedup_bam_spikein/{sample}_spikein.bam",
#         stat= "dedup_bam_spikein/{sample}_spikein.bam.stats.txt"
#     threads: 12
#     params:
#         spikein_chr = config["spike_in_chr"]
#     shell:
#         ## --threads flag failed
#         "(samtools view  -@ {threads} -hbS {input} {params.spikein_chr} | "
#         "samtools  sort  -@ {threads} -o {output.bam} && "
#         "samtools  index -@ {threads} {output.bam} && "
#         "samtools  stats -@ {threads} {output.bam} > {output.stat})"


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
        "(java -jar /cluster/tools/software/picard/2.10.9/picard.jar "
        "CollectInsertSizeMetrics M=0.01 I={input} O={output.txt} "
        "H={output.hist}) 2> {log}"

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
        "(java -jar /cluster/tools/software/picard/2.10.9/picard.jar "
        "CollectInsertSizeMetrics M=0.01 I={input} O={output.txt} "
        "H={output.hist}) 2> {log}"



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
        metrics_files=lambda wildcards: expand("fragment_size_virus/{sample}_insert_size_metrics.txt", sample=get_sample_ids_from_checkpoint(wildcards))
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


# ## spike-ins
# rule insert_size_spikein:
#     input:
#         "dedup_bam_spikein/{sample}_spikein.bam"
#     output:
#         txt = "fragment_size_spikein/{sample}_insert_size_metrics.txt",
#         hist = "fragment_size_spikein/{sample}_insert_size_histogram.pdf"
#     params:
#         pipeline_env = env_dir
#     log:
#         "logs/{sample}_picard_insert_size_spikein.log"
#     shell:
#         "(java -jar {params.pipeline_env}/share/picard-2.26.6-0/picard.jar "
#         "CollectInsertSizeMetrics M=0.05 I={input} O={output.txt} "
#         "H={output.hist}) 2> {log}"
