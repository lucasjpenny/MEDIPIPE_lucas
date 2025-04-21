#!/bin/bash
#SBATCH -p all              ## Specify SLURM partition for job submssion
#SBATCH -t 2-00:00:00
#SBATCH --mem=5G
#SBATCH -J submit_snakemake_%j
#SBATCH -o submit_snakemake_%j.out
#SBATCH -e submit_snakemake_%j.err

## submit the snakemake for sample list
## configure shell: full path to conda.sh
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate frag-pipeline

module load trim_galore
# module load samtools/1.14
module load samtools
module load picard/2.10.9
module load bwa
# module load bedtools/2.27.1
# module load R
module load fastqc
module load bowtie2

## now source conda and activate the right env
source ~/miniconda3/etc/profile.d/conda.sh 
conda activate frag-pipeline


## mkdir for cluster submission logs
## defined in .workflow/config/cluster_std_err.json
# cd  /path/to/test/Res
# mkdir -p logs_cluster

## unlock workdir just in case the folder locked accidently before
snakemake --snakefile /cluster/home/t116306uhn/workflows/MEDIPIPE_lucas/workflow/Snakefile  \
          --configfile /cluster/home/t116306uhn/workflows/MEDIPIPE_lucas/workflow/20250420_test_config_template.yaml \
          --unlock

## -p   partition to submit for SLURM
## --mem    request memory, specify as much as you can. As we tested, 60G can handle all real cfMeDIP-seq datasets so far.
## --jobs   ## number of samples in modified sample_tmplate.tsv files (independent steps will be scheduled per sample)
## sbatch -c :  maximal 12 threads per multithreading job by default, less -c INT  will be scaled down to INT


snakemake --snakefile /cluster/home/t116306uhn/workflows/MEDIPIPE_lucas/workflow/Snakefile \
          --configfile /cluster/home/t116306uhn/workflows/MEDIPIPE_lucas/workflow/20250420_test_config_template.yaml \
          --use-conda  --conda-prefix frag-pipeline \
          --cluster-config /cluster/home/t116306uhn/workflows/MEDIPIPE_lucas/workflow/config/cluster_std_err.json \
          --cluster "sbatch -p himem -c 12 --mem=60G -t 12:0:0 -o {cluster.std} -e {cluster.err}" \
          --latency-wait 60 --jobs 50 -p --rerun-incomplete #--dry-run --dag | dot -Tsvg > dag.svg

# conda deactivate