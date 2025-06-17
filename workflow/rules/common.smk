workdir: config['work_dir']      ## Configure Working Directory

import os
import pandas as pd

###########################################
## paths for pipeline and/or reference data
work_dir = config["work_dir"]
pipe_dir = config["pipe_dir"]
#env_dir = config["pipeline_env"]  #${CONDA_PREFIX} dosen't work
env_dir = os.getenv("CONDA_PREFIX")


##################################################q
## read in sample and corresponding fq files talbe
## and aggregatio talbe
SAMPLES = (
    pd.read_csv(config["samples"], sep="\t")
    .set_index("sample_id", drop=False)
    .sort_index()
)

## read in samples for aggregation
if config["aggregate"]:
    SAMPLES_AGGR = (
        pd.read_csv(config["samples_aggr"], sep="\t")
        .set_index("sample_id", drop=False)
        .sort_index()
    )
else:
    SAMPLES_AGGR = SAMPLES     ## must be defined


###############################
## read in refrence files' info
REF = pd.read_csv(config["ref_files"], sep="\t", header = None, index_col = 0)
blacklist = REF.loc["blacklist"][1]   ## ENCODE blacklist

def get_blacklist():
    return blacklist

def get_sample_ids_from_checkpoint(wildcards):
    # Read the CSV file produced by the checkpoint
    hpvviewer_df = checkpoints.extract_hpviewer_summary.get().output[0]
    df = pd.read_csv(hpvviewer_df)
    # Extract the sample IDs (assuming they are in the first column)
    sample_ids = df.iloc[:, 0].tolist()
    return sample_ids

#############################################
## get taget outputs based on the config file
## either for individual samples or aggregate
## all samples listed in sample_aggr.tsv !!!
#############################################

def get_rule_all_input():
    ## ensure extra env installed
    extra_env = "extra_env/all_extra_env_installed",

    ## fixed outputs
    #meth_qc = "aggregated/meth_qc.txt",
    aggr_qc = "aggregated/aggr_qc_report.html",
    meta_quant = "aggregated/meth_count.txt.gz",
    meth_filt = "autos_bfilt/meth_count_autos_bfilt.txt.gz",
    ####################
    # DELETE THIS
    fq_qc = expand("fastqc_pe/{samples}_R2_fastqc.zip", samples = SAMPLES["sample_id"]),
    meth_out = expand("meth_qc_quant/{samples}_count.txt", samples = SAMPLES["sample_id"]),
    frag_size = expand("fragment_size/{samples}_insert_size_metrics.txt", samples = SAMPLES["sample_id"]),
    fp_gc = expand("fragment_profile/{samples}_50_Granges.bed", samples = SAMPLES["sample_id"]),
    frag_size_secondary = expand("fragment_size_secondary/{samples}_insert_size_metrics.txt", samples = SAMPLES["sample_id"]),
    hpv_viewer_repeatmasker = expand("hpv_viewer_repeatmasker/{samples}/{samples}_HPV_profile.txt", samples = SAMPLES["sample_id"]),
    fusion_raw = expand("raw_bam_fusion/{samples}.bam", samples = SAMPLES["sample_id"]),
    frag_agg =  expand("fragment_size/fragment_length_summary.csv"),
    frag_agg_virus =  expand("fragment_size_virus/fragment_length_summary.csv"),
    end_motif_summary = expand("end_motif/MDS_scores.txt"),
    coverage = expand("HPV_coverage/{samples}_depths.txt", samples = SAMPLES["sample_id"]),
    coverage_filtered = expand("HPV_coverage/{samples}_above150bp_depths.txt",samples = SAMPLES["sample_id"]),
    binned_fragmentomics = expand("binned_bams_virus/{samples}_midpoints.bed",samples = SAMPLES["sample_id"]),

    return   frag_agg_virus + coverage + coverage_filtered + binned_fragmentomics#+ end_motif_summary + fq_qc + frag_size + frag_agg  # + hpv_viewer_repeatmasker
    ###################################
    ######################################
    ## aggregated outputs for SAMPLES_aggr
    ## paired-end and spike-in
    if config["aggregate"] and config["paired-end"] and config["spike_in"] and config["frag_profile"]:
        #mult_qc = "aggregated/QC_pe/multiqc_report.html",

        ## spike-ins
        spikein_mult_qc = "aggregated_spikein/QC_pe/multiqc_report.html",
        spikein_meth_qc = "aggregated_spikein/meth_qc.txt",
        spikein_meta_quant = "aggregated_spikein/meth_count.txt.gz",

        #fragment profiles
        fp_gc = "aggregated/fragment_profile_GC_corrected_1mb.tsv",        ## GC corrected fragment profile

        return  extra_env + aggr_qc + meta_quant + meth_filt + spikein_mult_qc + spikein_meth_qc + spikein_meta_quant + fp_gc

    ## paired-end and no spike-in
    elif config["aggregate"] and config["paired-end"] and config["spike_in"] == False and config["frag_profile"]:
        #mult_qc = "aggregated/QC_pe/multiqc_report.html",
        fp_gc = "aggregated/fragment_profile_GC_corrected_1mb.tsv",        ## GC corrected fragment profile

        return  extra_env + aggr_qc + meta_quant + meth_filt + fp_gc

    elif config["aggregate"] and config["paired-end"] and config["spike_in"] == False and config["frag_profile"] == False:
        #mult_qc = "aggregated/QC_pe/multiqc_report.html",
        return  extra_env + aggr_qc + meta_quant + meth_filt

    ## single-end
    elif config["aggregate"] and config["paired-end"] == False:
        #mult_qc = "aggregated/QC_se/multiqc_report.html",
        return extra_env + aggr_qc + meta_quant + meth_filt


    #####################################
    ## outputs for each individual sample
    ## paired-end and spike-in
    elif config["aggregate"] == False and config["paired-end"] and config["spike_in"] and config["frag_profile"]:
        fq_qc = expand("fastqc_pe/{samples}_R2_fastqc.zip", samples = SAMPLES["sample_id"]),
        meth_out = expand("meth_qc_quant/{samples}_count.txt", samples = SAMPLES["sample_id"]),
        meth_spikein = expand("meth_qc_quant_spikein/{samples}_count.txt", samples = SAMPLES["sample_id"]),
        frag_size = expand("fragment_size/{samples}_insert_size_metrics.txt", samples = SAMPLES["sample_id"]),
        fp_gc = expand("fragment_profile/{samples}_50_Granges.bed", samples = SAMPLES["sample_id"]),

        return extra_env + fq_qc + frag_size + meth_out + meth_spikein + fp_gc

    ## no frag_profile
    elif config["aggregate"] == False and config["paired-end"] and config["spike_in"] and config["frag_profile"] == False:
        fq_qc = expand("fastqc_pe/{samples}_R2_fastqc.zip", samples = SAMPLES["sample_id"]),
        meth_out = expand("meth_qc_quant/{samples}_count.txt", samples = SAMPLES["sample_id"]),
        meth_spikein = expand("meth_qc_quant_spikein/{samples}_count.txt", samples = SAMPLES["sample_id"]),
        frag_size = expand("fragment_size/{samples}_insert_size_metrics.txt", samples = SAMPLES["sample_id"]),

        return extra_env + fq_qc + frag_size + meth_out + meth_spikein

    ## paired-end without spike-in
    elif config["aggregate"] == False and config["paired-end"] and config["spike_in"] == False and config["frag_profile"]:
        fq_qc = expand("fastqc_pe/{samples}_R2_fastqc.zip", samples = SAMPLES["sample_id"]),
        meth_out = expand("meth_qc_quant/{samples}_count.txt", samples = SAMPLES["sample_id"]),
        frag_size = expand("fragment_size/{samples}_insert_size_metrics.txt", samples = SAMPLES["sample_id"]),
        fp_gc = expand("fragment_profile/{samples}_50_Granges.bed", samples = SAMPLES["sample_id"]),

        return extra_env + fq_qc + frag_size + meth_out + fp_gc

    ## paired-end without spike-in, and frag_profile
    elif config["aggregate"] == False and config["paired-end"] and config["spike_in"] == False and config["frag_profile"] == False:
        fq_qc = expand("fastqc_pe/{samples}_R2_fastqc.zip", samples = SAMPLES["sample_id"]),
        meth_out = expand("meth_qc_quant/{samples}_count.txt", samples = SAMPLES["sample_id"]),
        frag_size = expand("fragment_size/{samples}_insert_size_metrics.txt", samples = SAMPLES["sample_id"]),

        # I removed meth_out from this function!! and extra_env +
        return  fq_qc + frag_size 

    ## single-end
    elif config["aggregate"] == False and config["paired-end"] == False:
         fq_qc = expand("fastqc_se/{samples}_R2_fastqc.zip", samples = SAMPLES["sample_id"]),
         meth_out = expand("meth_qc_quant/{samples}_count.txt", samples = SAMPLES["sample_id"]),

         return extra_env + fq_qc + meth_out




############################
## other functions for input
#############################

###############################
##  get corresponding bwa_index
def get_bwa_index():
    if config["spike_in"]:
        #return REF.loc["bwa_idx_spikein"][1]
        return config["spike_idx"]
    else:
        return REF.loc["bwa_index"][1]

def get_bwa_fusion_index():
    if config["spike_in"]:
        #return REF.loc["bwa_idx_spikein"][1]
        return config["spike_idx"]
    else:
        return REF.loc["fusion_index"][1]

def get_bwa_shifted_index():
    if config["spike_in"]:
        #return REF.loc["bwa_idx_spikein"][1]
        return config["spike_idx"]
    else:
        return REF.loc["shifted"][1]


################################################################
##  get raw fastq files for rename to consistant wildcard.sample
def get_raw_fastq_se(wildcards):
    if config["paired-end"] == False:
        return SAMPLES.loc[wildcards.sample]["R1"].split(",")
    else:
        return ""

def get_raw_fastq_pe_R1(wildcards):
    if config["paired-end"]:
        return SAMPLES.loc[wildcards.sample]["R1"].split(",")
    else:
        return ""

def get_raw_fastq_pe_R2(wildcards):
    if config["paired-end"]:
        return SAMPLES.loc[wildcards.sample]["R2"].split(",")
    else:
        return ""



#######################################################
##  get renamed fastq for FASQC and barcode extraction
def get_renamed_fastq(wildcards):
    if config["paired-end"]:
        R1 = "renamed_fq/{}_R1.fastq.gz".format(wildcards.sample),
        R2 = "renamed_fq/{}_R2.fastq.gz".format(wildcards.sample),
        return R1 + R2
    else:
        return "renamed_fq/{}.fastq.gz".format(wildcards.sample)


#####################################################
## Get the HPV genotype (detected from Jingfeng's code)

def get_sample_hpv_genotype(sample):
    if config["paired-end"]:
        genotype = SAMPLES.loc[sample]["HPV_genotype"]
        return f"{config['hpv_dir']}/{genotype}.fasta"
    else:
        return ""


def get_sample_hpv_genotype_shifted(sample):
    """sample is just a plain Python string, e.g. "OPC_Pool-1_14_S14_L001"."""
    if config["paired-end"]:
        genotype = SAMPLES.loc[sample, "HPV_genotype"]
        return f"{config['hpv_dir']}/{genotype}L1toE1.fasta"
    else:
        return ""



def get_sample_hpv_info_shifted(sample):
    import os

    genotype = SAMPLES.loc[sample]["HPV_genotype"]
    fasta_path = os.path.join(config["hpv_dir"], f"{genotype}.fasta")

    with open(fasta_path, "r") as f:
        first_line = f.readline().strip()
        if first_line.startswith(">"):
            fasta_header = first_line[1:].split(" ")[0]  # remove '>' and take up to first space
        else:
            raise ValueError(f"FASTA file {fasta_path} doesn't start with a proper header line")

    region = f"{fasta_header}:3000-5000"
    return region



################################################
## get fastq for TRIM GALORE
## UMI extracted for paired-end reads, if exist
def get_fastq_4trim(wildcards):
    if config["paired-end"] and config["add_umi"]:
        R1 = "barcoded_fq_pe/{}_R1.fastq.gz".format(wildcards.sample),
        R2 = "barcoded_fq_pe/{}_R2.fastq.gz".format(wildcards.sample),
        return R1 + R2
    elif config["paired-end"] and config["add_umi"] == False:
        R1 = "renamed_fq/{}_R1.fastq.gz".format(wildcards.sample),
        R2 = "renamed_fq/{}_R2.fastq.gz".format(wildcards.sample),
        return R1 + R2
    elif config["paired-end"] == False and config["add_umi"]:
        return "barcoded_fq_se/{}.fastq.gz".format(wildcards.sample)
    else:
        return "renamed_fq/{}.fastq.gz".format(wildcards.sample)

##################################
## get trimmed fastq files for BWA
def get_trimmed_fastq(wildcards):
    if config["paired-end"]:
        R1 = "trimmed_fq/{}_R1_val_1.fq.gz".format(wildcards.sample),
        R2 = "trimmed_fq/{}_R2_val_2.fq.gz".format(wildcards.sample),
        return R1 + R2
    else:
        return "trimmed_fq/{}_trimmed.fq.gz".format(wildcards.sample)



########################################
## get dedup bam files for meth_qc_quant
def get_dedup_bam(wildcards):
    if config["paired-end"] and config["add_umi"]:
        return "dedup_bam_umi_pe/{}_dedup.bam".format(wildcards.sample)
    elif config["paired-end"] and config["add_umi"] == False:
        return "dedup_bam_pe/{}_dedup.bam".format(wildcards.sample)
    elif config["paired-end"] == False and config["add_umi"]:
        return "dedup_bam_umi_se/{}_dedup.bam".format(wildcards.sample)
    else:
        return "dedup_bam_se/{}_dedup.bam".format(wildcards.sample)

def get_dedup_bam_virus(wildcards):
    if config["paired-end"] and config["add_umi"]:
        return "dedup_bam_umi_pe_shifted/{}_dedup.bam".format(wildcards.sample)
    elif config["paired-end"] and config["add_umi"] == False:
        return "dedup_bam_umi_pe_shifted/{}_dedup.bam".format(wildcards.sample)
    elif config["paired-end"] == False and config["add_umi"]:
        return "dedup_bam_umi_pe_shifted/{}_dedup.bam".format(wildcards.sample)
    else:
        return "dedup_bam_umi_pe_shifted/{}_dedup.bam".format(wildcards.sample)


## spike-ins
def get_dedup_bam_spikein(wildcards):
    if config["paired-end"] and config["spike_in"]:
        return "dedup_bam_spikein/{}_spikein.bam".format(wildcards.sample)


#####################
## aggregaton #######
#####################

####################
## get FASTQC stats
def get_fastqc_stats():
    if config["paired-end"]:
        r1_raw  = expand("fastqc_pe/{samples}_R1_fastqc.zip", samples = SAMPLES_AGGR["sample_id"]),
        r2_raw  = expand("fastqc_pe/{samples}_R2_fastqc.zip", samples = SAMPLES_AGGR["sample_id"]),
        r1_trim = expand("fastqc_pe/{samples}_R1_val_1_fastqc.zip", samples = SAMPLES_AGGR["sample_id"]),
        r2_trim = expand("fastqc_pe/{samples}_R2_val_2_fastqc.zip", samples = SAMPLES_AGGR["sample_id"]),
        return r1_raw + r2_raw + r1_trim + r2_trim
    else:
         r1_raw  = expand("fastqc_se/{samples}_fastqc.zip", samples = SAMPLES_AGGR["sample_id"]),
         r1_trim = expand("fastqc_se/{samples}_trimmed_fastqc.zip", samples = SAMPLES_AGGR["sample_id"]),
         return r1_raw + r1_trim


######################
## get dedup bam stats
def get_dedup_bam_stats():
    if config["paired-end"] and config["add_umi"]:
        return expand("dedup_bam_umi_pe/{samples}_dedup.bam.stats.txt", samples = SAMPLES_AGGR["sample_id"])
    elif config["paired-end"] and config["add_umi"] == False:
        return expand("dedup_bam_pe/{samples}_dedup.bam.stats.txt", samples = SAMPLES_AGGR["sample_id"])
    elif config["paired-end"] == False and config["add_umi"]:
        return expand("dedup_bam_umi_se/{samples}_dedup.bam.stats.txt", samples = SAMPLES_AGGR["sample_id"])
    else:
        return expand("dedup_bam_se/{samples}_dedup.bam.stats.txt", samples = SAMPLES_AGGR["sample_id"])



############################
## get aggregated qc results
def get_aggr_qc_stats():
    if config["paired-end"]:
        fastqc = "aggregated/QC_pe/multiqc_data/multiqc_fastqc.txt",
        sam_stats = "aggregated/QC_pe/multiqc_data/multiqc_samtools_stats.txt",
        frag_stats = "aggregated/QC_pe/multiqc_data/multiqc_picard_insertSize.txt",
        meth_qc = "aggregated/meth_qc.txt",
        return fastqc + sam_stats + frag_stats + meth_qc
    else:
         fastqc = "aggregated/QC_se/multiqc_data/multiqc_fastqc.txt",
         sam_stats = "aggregated/QC_se/multiqc_data/multiqc_samtools_stats.txt",
         meth_qc = "aggregated/meth_qc.txt",
         return fastqc + sam_stats + meth_qc
