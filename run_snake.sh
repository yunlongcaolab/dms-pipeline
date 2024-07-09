#!/bin/bash

set -e

partitions=cpu1,cpu2,sugon,hygon
workdir=/lustre/grp/cyllab/share/DMS/$1

timestamp=$(date +%Y%m%d%H%M%S)

config=$workdir/config.yaml
snakefile=$workdir/Snakefile.$timestamp

slurm_script=$workdir/slurm_submit.$timestamp.sh

if [ -z $1 ]; then
    echo "Please provide the project name as input!"
    exit 1
fi

if [ ! -d $workdir ]; then
    echo "Workdir not found!"
    exit 1
fi

if [ ! -f $config ]; then
    echo "Config file not found!"
    exit 1
fi

mkdir -p $workdir/.old.files

# if the workdir contains old files, move them to .old.files

ls $workdir | grep ^Snakefile | xargs -I {} mv $workdir/{} $workdir/.old.files/
ls $workdir | grep ^slurm_submit | xargs -I {} mv $workdir/{} $workdir/.old.files/

cp Snakefile $snakefile

echo "Copied Snakefile: $snakefile"

# snakemake --snakefile $snakefile --configfile $config --cluster "$slurm_cmd" -j unlimited --default-resources cpu_per_task=1 mem_mb=2048 # for snakemake v7

slurm_cmd="sbatch --partition=$partitions -c {resources.cpu_per_task} --mem={resources.mem_mb} -J {rule}_{wildcards} --time=2400 -D $workdir -o {resources.stdout} -e {resources.stderr}" 

cat << EOF > $slurm_script
#!/bin/bash
#SBATCH -c 1
#SBATCH --mem=1g
#SBATCH --partition=$partitions
#SBATCH --time=2400
#SBATCH -o $workdir/slurm_submit_$timestamp.o.txt
#SBATCH -e $workdir/slurm_submit_$timestamp.e.txt
#SBATCH -D $workdir

python $(pwd)/scripts/sample_info_clean.py --config=$config --workdir=$workdir

snakemake --snakefile $snakefile --configfile $config --executor cluster-generic --cluster-generic-submit-cmd "$slurm_cmd" -j unlimited --default-resources cpu_per_task=1 mem_mb=2048
EOF

sbatch $slurm_script
echo "Submitted slurm job: $slurm_script"