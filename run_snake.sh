#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=1g
#SBATCH --partition=cpu1,cpu2,sugon,hygon

python scripts/sample_info_clean.py --config=config.yaml

# if scaat is available, use the following command (recommended by snakemake v8+). All slurm logs will be in .snakemake/slurm_logs
# snakemake --executor slurm -j unlimited --default-resources cpu_per_tasks=1 mem_mb=2048 time=1440 slurm_partition=cpu1,cpu2,hygon,sugon

# if scaat is banned, use generic executor. custom log files are used.
slurm_cmd="sbatch --partition=cpu1,cpu2,sugon,hygon -c {resources.cpu_per_task} --mem={resources.mem_mb} -J {rule}_{wildcards} --time=1440 -o {resources.stdout} -e {resources.stderr}" 
snakemake --executor cluster-generic --cluster-generic-submit-cmd "$slurm_cmd" -j unlimited --default-resources cpu_per_task=1 mem_mb=2048 # for snakemake v8+
# snakemake --cluster "$slurm_cmd" -j unlimited --default-resources cpu_per_task=1 mem_mb=2048 # for snakemake v7
