#!/bin/bash

set -e

partitions=cpu1,cpu2,sugon,hygon

timestamp=$(date +%Y%m%d%H%M%S)

workdir=$(pwd)
rule=$1
config=$workdir/config.yaml
snakefile=$workdir/Snakefile

slurm_script=$workdir/slurm_test/slurm_submit.$timestamp.sh

if [ ! -f $config ]; then
    echo "Config file not found!"
    exit 1
fi

slurm_cmd="sbatch --partition=$partitions -c {resources.cpu_per_task} --mem={resources.mem_mb} -J {rule}_{wildcards} --time=2400 -D $workdir -o {resources.stdout} -e {resources.stderr}" 

cat << EOF > $slurm_script
#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=1g
#SBATCH --partition=$partitions
#SBATCH --time=2400
#SBATCH -o $workdir/slurm_test/slurm_submit_%j_$timestamp.o.txt
#SBATCH -e $workdir/slurm_test/slurm_submit_%j_$timestamp.e.txt
#SBATCH -D $workdir

python $(pwd)/scripts/sample_info_clean.py --config=$config --workdir=$workdir

snakemake $rule --snakefile $snakefile --configfile $config --executor cluster-generic --cluster-generic-submit-cmd "$slurm_cmd" -j unlimited --default-resources cpu_per_task=1 mem_mb=2048
EOF

sbatch $slurm_script
echo "Submitted slurm job: $slurm_script"
