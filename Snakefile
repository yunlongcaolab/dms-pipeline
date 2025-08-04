import os, glob
import pandas as pd
import yaml

from pathlib import Path

with open(os.path.join(config['output'], "_tasks.yaml"), "r") as f:
    TASKS = yaml.safe_load(f)
    if config['batches'] != 'all':
        print(f"Using specified batches {config['batches']}", flush=True)
        BATCHES = config['batches']
    else:
        print("Using all batches", flush=True)
        BATCHES = set()
        for task in TASKS['escape_calc']:
            BATCHES.add(task)
        for sample in TASKS['barcode_count']:
            BATCHES.add(TASKS['barcode_count'][sample]['batch'])

if 'default_libinfo' in config:
    for key in config['libinfo']:
        for k, v in config['default_libinfo'].items():
            if k not in config['libinfo'][key]:
                config['libinfo'][key][k] = v

rule all:
    input:
        os.path.join(config['output'], "barcode_count_stat.csv"),
        expand(os.path.join(config['output'], "escape_summary/{batch}/QCstat.pdf"), batch=TASKS['escape_calc'].keys()),
        expand(os.path.join(config['output'], "barcode_count/{sample}/overrepresented_sequences.txt"), sample=TASKS["barcode_count"].keys()),
        [os.path.join(config['output'], f"sort_seq/{batch}/{library}/variant_scores.csv") for batch in TASKS['sort_seq'] for library in TASKS['sort_seq'][batch]]

targets = os.listdir(config['raw_lib'])
tables = []

for tg in targets:
    libs = os.listdir(os.path.join(config['raw_lib'], tg))
    for lib in libs:
        tables.append(os.path.join(config['output'], 'library_tables', tg, lib, 'variant_table.csv'))

if 'library_merge' in config:
    for tg in config['library_merge']:
        for lib in config['library_merge'][tg]:
            tables.append(os.path.join(config['output'], 'library_tables', tg, lib, 'variant_table.csv'))

rule all_tables:
    input:
        tables

rule generate_target_ref:
    input:
        plasmid = os.path.join(config['target_ref'], '{target}/plasmid.txt'),
        template = lambda wc: os.path.join(config['target_ref'], wc.target, 'template.txt') if os.path.exists(os.path.join(config['target_ref'], wc.target, 'template.txt')) else config['default_template'],
    output:
        minimap_ref = os.path.join(config['target_ref'], '{target}/{target}.fasta'),
        gb = os.path.join(config['target_ref'], '{target}/pacbio_amplicons.gb'),
        specs = os.path.join(config['target_ref'], '{target}/parse_specs.yaml'),
        wt_seq = os.path.join(config['wt_seqs'], '{target}.fasta') if 'wt_seqs' in config else os.path.join(config['target_ref'], '..','wt_seqs/{target}.fasta')
    resources:
        stdout = lambda wc: os.path.join(config['output'], f"logs/generate_target_ref/{wc.target}_stdout.txt"),
        stderr = lambda wc: os.path.join(config['output'], f"logs/generate_target_ref/{wc.target}_stderr.txt")
    script:
        f"{config['pipeline']}/scripts/generate_target_ref.py"

rule ccs_align:
    input:
        reads = lambda wc: (Path(config['raw_lib']) / wc.target / wc.library).rglob('*.f*q.gz'),
        ref = os.path.join(config['target_ref'], '{target}/{target}.fasta')
    output:
        os.path.join(config['output'], 'library_tables/{target}/{library}/alignments.bam')
    resources:
        stdout = lambda wc: os.path.join(config['output'], f"logs/ccs_align/{wc.target}/{wc.library}_stdout.txt"),
        stderr = lambda wc: os.path.join(config['output'], f"logs/ccs_align/{wc.target}/{wc.library}_stderr.txt"),
        cpu_per_task = config['cpu_per_task']
    shell:
        f"mkdir -p {os.path.join(config['output'], 'library_tables/{wildcards.target}/{wildcards.library}')} && "
        f"minimap2 -a -A5 -B7 -O16 -E2 --end-bonus=23 --secondary=no --cs -t {config['cpu_per_task']} {{input.ref}} {{input.reads}} | samtools view -bS > {{output}}"

rule library_table:
    input:
        aln = os.path.join(config['output'], 'library_tables/{target}/{library}/alignments.bam'),
        ref = os.path.join(config['target_ref'], '{target}/pacbio_amplicons.gb'),
        specs = os.path.join(config['target_ref'], '{target}/parse_specs.yaml'),
        reads = lambda wc: (Path(config['raw_lib']) / wc.target / wc.library).rglob('*.f*q.gz')
    output:
        os.path.join(config['output'], 'library_tables/{target}/{library}/variant_table.csv'),
        os.path.join(config['output'], 'library_tables/{target}/{library}/consensus.csv.gz'),
        os.path.join(config['output'], 'library_tables/{target}/{library}/processed_ccs.csv.gz'),
        os.path.join(config['output'], 'library_tables/{target}/{library}/filtered_ccs.csv.gz'),
    resources:
        stdout = lambda wc: os.path.join(config['output'], f"logs/library_table/{wc.target}/{wc.library}_stdout.txt"),
        stderr = lambda wc: os.path.join(config['output'], f"logs/library_table/{wc.target}/{wc.library}_stderr.txt")
    script:
        f"{config['pipeline']}/scripts/parse_ccs_alignments.py"

if 'library_merge' in config:
    ruleorder: library_merge > library_table
    rule library_merge:
        input:
            tables=lambda wc: [
                os.path.join(config['output'], f'library_tables/{wc.target}/{lib}/variant_table.csv') for lib in config['library_merge'][wc.target][wc.library]
            ],
            ref = lambda wc: (
                os.path.join(config['wt_seqs'], wc.target + '.fasta') 
                if 'wt_seqs' in config else 
                    os.path.join(config['target_ref'], '..', 'wt_seqs', wc.target + '.fasta')
            )
        output:
            os.path.join(config['output'], 'library_tables/{target}/{library}/variant_table.csv')
        resources:
            stdout = lambda wc: os.path.join(config['output'], f"logs/library_merge/{wc.target}/{wc.library}_stdout.txt"),
            stderr = lambda wc: os.path.join(config['output'], f"logs/library_merge/{wc.target}/{wc.library}_stderr.txt")
        wildcard_constraints:
            target='|'.join(config['library_merge'].keys()),
            library='|'.join([lib for lib_list in config['library_merge'].values() for lib in lib_list])
        script:
            f"{config['pipeline']}/scripts/merge_library_table.py"

# rule library_QC_plot:
#     input:
#         os.path.join(config['output'], 'library_tables/{target}/{library}_variant_table.csv'),
#         os.path.join(config['output'], 'library_tables/{target}/{library}_stats.csv')
#     output:
#         os.path.join(config['output'], 'library_tables/{target}/{library}_QCplot.pdf')
#     resources:
#         stdout = os.path.join(config['output'], "logs/library_QC_plot/{target}/{library}_stdout.txt"),
#         stderr = os.path.join(config['output'], "logs/library_QC_plot/{target}/{library}_stderr.txt")
#     script:
#         f"{config['pipeline']}/scripts/library_QC_plot.py"

rule barcode_count:
    input:
        reads = lambda wildcards: TASKS['barcode_count'][wildcards.sample]["fastq_files"],
        table = lambda wildcards: os.path.join(config['output'], 'library_tables', config['libinfo'][TASKS['barcode_count'][wildcards.sample]['library']]['target'], TASKS['barcode_count'][wildcards.sample]['library'], 'variant_table.csv'),
        script = os.path.join(config['pipeline'], 'scripts', 'barcode_count.py')
    output:
        os.path.join(config['output'], "barcode_count/{sample}/counts.csv"),
        os.path.join(config['output'], "barcode_count/{sample}/barcode_count_info.yaml")
    params:
        library = lambda wildcards: TASKS['barcode_count'][wildcards.sample]['library'],
        batch = lambda wildcards: TASKS['barcode_count'][wildcards.sample]['batch'],
        fastq_files = lambda wildcards: ' '.join(TASKS['barcode_count'][wildcards.sample]["fastq_files"]),
        bclen = lambda wildcards: config['libinfo'][TASKS['barcode_count'][wildcards.sample]['library']]['bclen'],
        # table = lambda wildcards: os.path.join(config['output'], 'library_tables', config['libinfo'][TASKS['barcode_count'][wildcards.sample]['library']]['target'], TASKS['barcode_count'][wildcards.sample]['library']+'_variant_table.csv'),
    resources:
        stdout = lambda wc: os.path.join(config['output'], f"logs/barcode_count/{wc.sample}/stdout.txt"),
        stderr = lambda wc: os.path.join(config['output'], f"logs/barcode_count/{wc.sample}/stderr.txt")
    wildcard_constraints:
        sample='|'.join(TASKS["barcode_count"].keys())
    shell:
        f'''python {{input.script}} --sample={{wildcards.sample}} --batch={{params.batch}} --library={{params.library}} \
                                       -i {{params.fastq_files}} \
                                       -o {os.path.join(config['output'], 'barcode_count/{wildcards.sample}')} \
                                       -b {{params.bclen}} \
                                       -t {{input.table}} \
                                       --allowed_lowq={config['barcode_count']['allowed_lowq']} \
                                       --lowq={config['barcode_count']['lowq']} \
                                       --max_distance={config['barcode_count']['max_distance']} \
                                       --min_dist_diff={config['barcode_count']['min_dist_diff']}'''

def get_ref_merge_input(wildcards):
    files = []
    for sample in TASKS['ref_merge'][wildcards.target_ref]:
        for file in ["counts.csv", "barcode_count_info.yaml"]:
            files.append(os.path.join(config['output'], "barcode_count", sample, file))
    return files

rule ref_merge:
    input: get_ref_merge_input
    output:
        os.path.join(config['output'], "barcode_count/{target_ref}/counts.csv"),
        os.path.join(config['output'], "barcode_count/{target_ref}/barcode_count_info.yaml")
    resources:
        stdout = lambda wc: os.path.join(config["output"], f"logs/ref_merge/{wc.target_ref}/stdout.txt"),
        stderr = lambda wc: os.path.join(config["output"], f"logs/ref_merge/{wc.target_ref}/stderr.txt")
    wildcard_constraints:
        target_ref='|'.join(TASKS["ref_merge"].keys())
    script: 
        f"{config['pipeline']}/scripts/ref_merge.py"

rule barcode_count_stat:
    input:
        expand(os.path.join(config['output'], "barcode_count/{sample}/barcode_count_info.yaml"), sample=TASKS['barcode_count'].keys()),
        expand(os.path.join(config['output'], "barcode_count/{sample}/barcode_count_info.yaml"), sample=TASKS['ref_merge'].keys())
    output:
        os.path.join(config['output'], "barcode_count_stat.csv")

    resources:
        stdout = os.path.join(config["output"], "logs/barcode_count_stat/stdout.txt"),
        stderr = os.path.join(config["output"], "logs/barcode_count_stat/stderr.txt")
    script:
        f"{config['pipeline']}/scripts/barcode_count_stat.py"

rule fastqc:
    input:
        reads = lambda wildcards: TASKS['barcode_count'][wildcards.sample]["fastq_files"]
    output:
        os.path.join(config['output'], "barcode_count/{sample}/overrepresented_sequences.txt")
    resources:
        stdout = lambda wc: os.path.join(config['output'], f"logs/fastqc/{wc.sample}/stdout.txt"),
        stderr = lambda wc: os.path.join(config['output'], f"logs/fastqc/{wc.sample}/stderr.txt")
    params:
        fastqc_dir = os.path.join(config['output'], "barcode_count/{sample}/fastqc"),
    wildcard_constraints:
        sample='|'.join(TASKS["barcode_count"].keys())
    shell:
        """
        shopt -s extglob
        mkdir -p {params.fastqc_dir} && \
        fastqc --extract -o {params.fastqc_dir} {input.reads} && \
        rm -rf {params.fastqc_dir}/*@(_2|_R2|.R2)_fastqc && \
        awk '/>>Overrepresented sequences/,/>>END_MODULE/' \
            {params.fastqc_dir}/*_fastqc/fastqc_data.txt | \
            grep -v ">>" > {output} || true
        """

rule escape_calc:
    input:
        lambda wildcards: [
            os.path.join(config['output'], f"barcode_count/{wildcards.sample}/counts.csv"),
            os.path.join(config['output'], f"barcode_count/{TASKS['escape_calc'][wildcards.batch][wildcards.sample]['ref']}/counts.csv"),
            os.path.join(config['output'], f"barcode_count/{wildcards.sample}/barcode_count_info.yaml"),
            os.path.join(config['output'], f"barcode_count/{TASKS['escape_calc'][wildcards.batch][wildcards.sample]['ref']}/barcode_count_info.yaml"),
        ]
    params:
        library = lambda wildcards: TASKS['escape_calc'][wildcards.batch][wildcards.sample]['library'],
        antibody = lambda wildcards: TASKS['escape_calc'][wildcards.batch][wildcards.sample]['antibody'],
        ref = lambda wildcards: TASKS['escape_calc'][wildcards.batch][wildcards.sample]['ref'],
        table = lambda wildcards: os.path.join(config['output'], 'library_tables', config['libinfo'][TASKS['escape_calc'][wildcards.batch][wildcards.sample]['library']]['target'], TASKS['escape_calc'][wildcards.batch][wildcards.sample]['library'], 'variant_table.csv'),
        wt_seq = lambda wildcards: (
            os.path.join(config['wt_seqs'], config['libinfo'][TASKS['escape_calc'][wildcards.batch][wildcards.sample]['library']]['target'] + '.fasta') 
            if 'wt_seqs' in config else 
                os.path.join(config['target_ref'], '..', 'wt_seqs', 
                config['libinfo'][TASKS['escape_calc'][wildcards.batch][wildcards.sample]['library']]['target'] + '.fasta')
        )

    output:
        os.path.join(config['output'], "escape_calc/{batch}/{sample}/single_mut_escape_scores.csv"),
        os.path.join(config['output'], "escape_calc/{batch}/{sample}/site_escape_scores.csv"),
        os.path.join(config['output'], "escape_calc/{batch}/{sample}/variant_escape_scores.csv"),
        os.path.join(config['output'], "escape_calc/{batch}/{sample}/calc_escape_info.yaml"),
        os.path.join(config['output'], "escape_calc/{batch}/{sample}/calc_escape_stat.yaml")
    resources:
        stdout = lambda wc: os.path.join(config['output'], f"logs/escape_calc/{wc.batch}/{wc.sample}/stdout.txt"),
        stderr = lambda wc: os.path.join(config['output'], f"logs/escape_calc/{wc.batch}/{wc.sample}/stderr.txt")
    script:
        f"{config['pipeline']}/scripts/calc_escape_score.py"

rule escape_batch_summary:
    input:
        lambda wildcards: [
            *[os.path.join(config['output'], "escape_calc", wildcards.batch, sample, "site_escape_scores.csv") for sample in TASKS['escape_calc'][wildcards.batch]],
            *[os.path.join(config['output'], "escape_calc", wildcards.batch, sample, "calc_escape_stat.yaml") for sample in TASKS['escape_calc'][wildcards.batch]],
        ]
    output:
        expand(os.path.join(config['output'], "escape_summary/{{batch}}/{file}"), file=['stat.csv', 'site_escape_scores.csv'])
    resources:
        stdout = lambda wc: os.path.join(config["output"], f"logs/escape_summary/{wc.batch}/summary_stdout.txt"),
        stderr = lambda wc: os.path.join(config["output"], f"logs/escape_summary/{wc.batch}/summary_stderr.txt")
    script:
        f"{config['pipeline']}/scripts/escape_batch_summary.py"

rule escape_batch_plot:
    input:
        os.path.join(config['output'], "escape_summary/{batch}/site_escape_scores.csv"),
        os.path.join(config['output'], "escape_summary/{batch}/stat.csv"),
        os.path.join(config['output'], "barcode_count_stat.csv")
    output:
        expand(os.path.join(config['output'], "escape_summary/{{batch}}/{status}/site_escape_{agg}_{model}."+('png' if config.get('save_png', False) else 'pdf')), agg=['total', 'mean'], model=['model', 'single'], status=['pass', 'fail']),
        os.path.join(config['output'], "escape_summary/{batch}/QCstat.pdf")
    resources:
        stdout = lambda wc: os.path.join(config["output"], f"logs/escape_summary/{wc.batch}/plot_stdout.txt"),
        stderr = lambda wc: os.path.join(config["output"], f"logs/escape_summary/{wc.batch}/plot_stderr.txt")
    script:
        f"{config['pipeline']}/scripts/escape_batch_plot.py"

rule all_sort_seq:
    input:
        [os.path.join(config['output'], f"sort_seq/{batch}/{library}/variant_scores.csv") for batch in TASKS['sort_seq'] for library in TASKS['sort_seq'][batch]]
        # [os.path.join(config['output'], f"sort_seq/{batch}/{library}/single_mut_scores.csv") for batch in TASKS['sort_seq'] for library in TASKS['sort_seq'][batch]]

rule sort_seq:
    input:
        lambda wc: ([
            os.path.join(config['output'], f'barcode_count/{bin_info["sample"]}/counts.csv') for bin_info in TASKS['sort_seq'][wc.batch][wc.library].values()
        ])
    output:
        variant = os.path.join(config['output'], 'sort_seq/{batch}/{library}/variant_scores.csv'),
    params:
        bins = lambda wc: TASKS['sort_seq'][wc.batch][wc.library],
        table = lambda wc: os.path.join(config['output'], 'library_tables', config['libinfo'][wc.library]['target'], wc.library, 'variant_table.csv')
    resources:
        stdout = lambda wc: os.path.join(config['output'], f"logs/sort_seq/{wc.batch}/{wc.library}/stdout.txt"),
        stderr = lambda wc: os.path.join(config['output'], f"logs/sort_seq/{wc.batch}/{wc.library}/stderr.txt"),
        cpu_per_task = config['cpu_per_task']
    script:
        f"{config['pipeline']}/scripts/sort_seq.py"

