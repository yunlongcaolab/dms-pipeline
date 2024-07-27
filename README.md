
<h1><p align='center'>dms-pipeline: a Pipeline for High-throughput Antibody DMS Data Processing</p></h1>

<p align='center'>
  English | <a href='./README.zh-CN.md'>简体中文</a> 
</p>

This pipeline is designed for the processing and analysis of high-throughput deep mutational scanning (DMS) data. The upstream part of the analysis follows the idea of [pipeline from J. Bloom Lab](https://github.com/jbloomlab/SARS-CoV-2-RBD_DMS) with some custom modification to fit the high-throughput data.

You can refer to the published articles [here](#Citation) from our lab for more information. 

Notably, because some of the data presented in the articles may not be processed by the latest pipeline, the detailed results shown in the articles (including their corresponding Github repo) may not be identical with the results generated by this optimized pipeline. However, there should not be significant differences. 

Using the latest version of the pipeline in this repo is recommended.

## Installation
### Prepare conda

Most part of this pipeline is written in Python, and we will use conda to setup the environment. For accelerated env solving, using [Miniforge](https://github.com/conda-forge/miniforge) is highly recommended. If so, use `mamba` instead of `conda` in the following commands.

First, make sure your conda (or mamba) is ready by checking available environments.

```bash
conda env list
```

### Create environment

We provide an `environment.yml` file within this repo. Please clone this repo to your work directory and run

```bash
conda env create -f environment.yml
```

This will build a environment named `dms-pipeline`. Then, activate it

```bash
conda activate dms-pipeline
```

### Compile C/C++ code (optional)

A minor part of this pipeline is written in C/C++ for performance. These code provides optional features which are not necessary for basic analyses. If you need them, enter the `scripts` directory and compile them. Currently, the C/C++ code is very simple so there should not be complex issues on dependencies.

You will need to install `pybind11` by `pip` or `conda` to export the binaries as python modules.

```bash
conda install -c conda-forge pybind11

cd scripts
chmod +x compile.sh && ./compile.sh
```

If you failed in compiling these code, try installing a newer version of GCC using conda:

```bash
conda install -c conda-forge compilers
```

## Usage
### Configuration

See the example `config.yaml`. You should replace the items in it with your custom paths and other necessary information.
Note the following parts:

#### Relevant directories

#### libinfo

You need to specify the name and antigen (target) name of each library. The library name should be consistent with the directory name in `library_PacBio`.

All items specified in "default_libinfo" will be used as default values for all libraries. You can override them in the "libinfo" section.

#### Configuration of the parameters of the pipeline

### Files


We recommend the following directory structure to use this pipeline. First, build a directory named 'SomeProj' in any place you like, and clone this repository at another place (split code and data). 

Then, make the structure as follows:

```
SomeProj
├── NGS_raw
│   ├── batch1
│   ├── batch2
│   └── ...
├── NGS_sample_info
│   ├── batch1
│   │   └── sample_info.csv
│   ├── batch2
│   │   └── sample_info.csv
│   └── ...
├── PacBio_raw (raw PacBio data)
│   ├── batch1
│   ├── batch2
│   └── ...
|── library_PacBio (usually link from PacBio_raw)
|   ├── antigen1
|   │   ├── library1
|   │   │   ├── library1_run1.fq.gz
|   │   │   └── ...    
|   │   └── library2
|   └── ...
├── reference
│   ├── wt_seqs
│   │   ├── antigen1.fasta
│   │   ├── antigen2.fasta
│   │   └── ...
│   ├── filters
│   │   ├── antigen1
│   │   │   ├── filter1.csv
│   │   │   └── ...
│   │   └── ...
│   ├── target_ref
│   │   ├── antigen1
│   │   │   ├── plasmid.txt
│   │   │   └── template.txt
│   │   └── ...
│   └── ...
├── processed
│   └── ... (to be generated by the pipeline)
└── config.yaml (make sure you have written the correct paths in the config file)

dms-pipeline (this repo)
├── scripts
│   └── ...
├── run_snake.sh
├── Snakefile
└── ...
```

### Run the pipeline

Check the `run_snake.sh` script and run it to start the pipeline.

```bash
./run_snake.sh <directory> <rule>
```

Use `all` rule to process all data according to the configuration; use `all_tables` to run the PacBio part only.

It's OK to run the pipeline with other directory structures, as long as you write the `config.yaml` according.

### Downstream analyses

If you collected enough antibody DMS profiles for an antigen or several closely related antigens, you can conduct downstream analyses for it to get a comprehensive understand of its BCR epitopes.

We put the core code for the downstream analyses as another repository as a package named [HADI](https://github.com/yunlongcaolab/hadi), short for <ins>**H**</ins>igh-throughput <ins>**A**</ins>ntibody <ins>**D**</ins>MS <ins>**I**</ins>ntegration.

The script for downstream analyse using HADI is not integrated into the Snakemake pipeline because it is not always be conducted. 

You should properly write HADI config YAML files and run the `scripts/dms_integrate.py` script **(TO COMPLETE)**.

## TODO

- Include scripts for library table construction.
- Include Sort-seq and Tite-seq.
- Integrate Quality Control (QC) for DMS experiments in the pipeline, instead of HADI.

## Citation

Please cite the following papers if this pipeline is helpful:
- [Cao et al. Nature 2022](https://doi.org/10.1038/s41586-022-04980-y)
- [Cao et al. Nature 2023](https://doi.org/10.1038/s41586-022-05644-7)
- [Yisimayi et al. Nature 2024](https://doi.org/10.1038/s41586-023-06753-7)

You can find more information on [the official website of Cao lab](https://yunlongcaolab.com).
