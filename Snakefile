import os, glob
import pandas as pd
import yaml

with open(os.path.join(config['output'], "_tasks.yaml"), "r") as f:
    TASKS = yaml.safe_load(f)
    if config['batches'] != 'all':
        print(f"Using specified batches {config['batches']}", flush=True)
        BATCHES = config['batches']
    else:
        print("Using all batches", flush=True)
        BATCHES = list(TASKS['escape_calc'].keys())

rule all:
    input:
        expand(os.path.join(config['output'], "escape_summary/{batch}/QCstat.pdf"), batch=BATCHES),

rule barcode_count:
    input:
        lambda wildcards: TASKS['barcode_count'][wildcards.sample]["fastq_files"]
    output:
        os.path.join(config['output'], "barcode_count/{sample}/counts.csv"),
        os.path.join(config['output'], "barcode_count/{sample}/barcode_count_info.yaml")
    params:
        library = lambda wildcards: TASKS['barcode_count'][wildcards.sample]['library'],
        batch = lambda wildcards: TASKS['barcode_count'][wildcards.sample]['batch'],
        fastq_files = lambda wildcards: ' '.join(TASKS['barcode_count'][wildcards.sample]["fastq_files"]),
        bclen = lambda wildcards: config['libinfo'][TASKS['barcode_count'][wildcards.sample]['library']]['bclen'],
        table = lambda wildcards: os.path.join(config['output'], 'library_tables', config['libinfo'][TASKS['barcode_count'][wildcards.sample]['library']]['target'], TASKS['barcode_count'][wildcards.sample]['library']+'_variant_table.csv'),
    resources:
        stdout = lambda wc: os.path.join(config['output'], f"logs/barcode_count/{wc.sample}/stdout.txt"),
        stderr = lambda wc: os.path.join(config['output'], f"logs/barcode_count/{wc.sample}/stderr.txt")
    wildcard_constraints:
        sample='|'.join(TASKS["barcode_count"].keys())
    shell:
        f'''python {config['pipeline']}/scripts/barcode_count.py --sample={{wildcards.sample}} --batch={{params.batch}} --library={{params.library}} \
                                       -i {{params.fastq_files}} \
                                       -o {os.path.join(config['output'], 'barcode_count/{wildcards.sample}')} \
                                       -b {{params.bclen}} \
                                       -t {{params.table}} \
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
        expand(os.path.join(config['output'], "escape_summary/{{batch}}/{status}/site_escape_{agg}_{model}.pdf"), agg=['total', 'mean'], model=['model', 'single'], status=['pass', 'fail']),
        expand(os.path.join(config['output'], "escape_summary/{{batch}}/QCstat.pdf"))
    resources:
        stdout = lambda wc: os.path.join(config["output"], f"logs/escape_summary/{wc.batch}/plot_stdout.txt"),
        stderr = lambda wc: os.path.join(config["output"], f"logs/escape_summary/{wc.batch}/plot_stderr.txt")
    script:
        f"{config['pipeline']}/scripts/escape_batch_plot.py"