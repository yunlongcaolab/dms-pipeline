#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=1g
#SBATCH --partition=cpu1,cpu2,sugon,hygon

# check config file exists
if [ ! -f $1 ]; then
    echo "Snakefile not found!"
    exit 1
fi

snakefile=$1
workdir=$(dirname $snakefile)
config=$workdir/config.yaml

python scripts/sample_info_clean.py --config=$config --workdir=$workdir

slurm_cmd="sbatch --partition=cpu1,cpu2,sugon,hygon -c {resources.cpu_per_task} --mem={resources.mem_mb} -J {rule}_{wildcards} --time=1440 -D $workdir -o {resources.stdout} -e {resources.stderr}" 
snakemake --snakefile $snakefile --directory $workdir --executor cluster-generic --cluster-generic-submit-cmd "$slurm_cmd" -j unlimited --default-resources cpu_per_task=1 mem_mb=2048 # for snakemake v8+
# snakemake --snakefile $snakefile --directory $workdir --cluster "$slurm_cmd" -j unlimited --default-resources cpu_per_task=1 mem_mb=2048 # for snakemake v7
