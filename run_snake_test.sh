#!/bin/bash

set -e

partitions=cpu1,cpu2,cpu3,hygon

timestamp=$(date +%Y%m%d%H%M%S)

pipeline=$(pwd)
workdir=$(pwd)
rule=$1
config=$workdir/config.yaml
snakefile=$workdir/Snakefile

num_proc_tasks=8

slurm_script_tasks=$workdir/slurm_test/slurm_tasks.$timestamp.sh
slurm_script=$workdir/slurm_test/slurm_submit.$timestamp.sh

mkdir -p $workdir/slurm_test

if [ ! -f $config ]; then
    echo "Config file not found!"
    exit 1
fi

cat << EEOF > $slurm_script_tasks
#!/bin/bash
#SBATCH -c $num_proc_tasks
#SBATCH --mem=32g
#SBATCH --partition=$partitions
#SBATCH -o $workdir/slurm_test/slurm_tasks_%j_$timestamp.o.txt
#SBATCH -e $workdir/slurm_test/slurm_tasks_%j_$timestamp.e.txt
#SBATCH -D $workdir

set -e
python $pipeline/scripts/sample_info_clean.py --config=$config --workdir=$workdir

slurm_cmd="sbatch --partition=$partitions -c {resources.cpu_per_task} --mem={resources.mem_mb} -J {rule}_{wildcards} --time=2400 -D $workdir -o {resources.stdout} -e {resources.stderr}" 

cat << EOF > $slurm_script
#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=1g
#SBATCH --partition=$partitions
#SBATCH -o $workdir/slurm_test/slurm_submit_%j_$timestamp.o.txt
#SBATCH -e $workdir/slurm_test/slurm_submit_%j_$timestamp.e.txt
#SBATCH -D $workdir


snakemake $rule --snakefile $snakefile --config pipeline=$pipeline --configfile $config --executor cluster-generic --cluster-generic-submit-cmd "\$slurm_cmd" -j unlimited --default-resources cpu_per_task=1 mem_mb=2048
EOF

sbatch $slurm_script
echo "Submitted slurm job: $slurm_script"
EEOF

sbatch $slurm_script_tasks