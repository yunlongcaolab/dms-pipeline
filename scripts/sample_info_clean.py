import pandas as pd
import sys, os, re
import yaml
import argparse
import time

from pathlib import Path
import glob
from sort_seq import load_expression_info

from concurrent.futures import ProcessPoolExecutor
import logging
import fnmatch, re

def get_fq_files(all_candidate_files, batch, sample, suffix, R2_pattern, is_cleaned):
    logger = logging.getLogger(__name__)

    patterns = [
        f"{sample}."+suffix,
        f"{sample}_"+suffix,
        f"{sample}"+suffix,
        f"*{sample}."+suffix,
        f"*{sample}_"+suffix,
        f"*{sample}"+suffix,
    ] # search pattern from all_candidate_files until one file is found
    for pattern in patterns:
        regex = re.compile(fnmatch.translate(pattern))
        fq_files = [x for x in all_candidate_files if regex.match(Path(x).name)]
        if len(fq_files) > 0:
            break

    if is_cleaned:
        fq_files = [x for x in fq_files if not '/CleanData/' in x]

    if len(fq_files) == 0:
        logger.warning(
            f"Warning: {batch} - {sample}. FASTQ not found for pattern '*{sample}{suffix}'. Try alternative path."
        )
        alt_pattern = re.compile(
            fnmatch.translate(os.path.join(f"*{sample}*", suffix if suffix[0] == "*" else "*" + suffix))
        )
        fq_files = [x for x in all_candidate_files if alt_pattern.match(x)]

    # if QC exist, remove all non-QC files
    for x in fq_files:
        if '/QC/' in x:
            fq_files = [x for x in fq_files if '/QC/' in x]
            break

    for x in fq_files:
        if '/Cleandata/' in x:
            fq_files = [x for x in fq_files if '/Cleandata/' in x]
            break

    if len(fq_files) == 0:
        logger.warning(
            f"Warning: {batch} - {sample}. FASTQ not found in alternative path. Skipped."
        )
        return ""
    elif "," in "".join(fq_files):
        logger.error(
            f"Error: {batch} - {sample}. Comma in fastq file path."
        )
        raise RuntimeError(f"Comma in fastq file path: {batch} - {sample}")
    elif len(fq_files) > 1:
        logger.warning(
            f"Warning: {batch} - {sample}. Multiple FASTQ found. {','.join(fq_files)}. Try removing R2 or redundant files."
        )
        fq_files = [x for x in fq_files if not R2_pattern.search(x)]

        logger.info(f"    {len(fq_files)} files found after removing. {','.join(fq_files)}.")
        return ",".join(fq_files)
    else:
        return fq_files[0]

def parse_batch(sample_info, raw, batch, necessary_columns, sample_name_replace, exclude_antibody_names, libraries, suffix, R2_pattern):
    logger = logging.getLogger(__name__)
    csvfile = sample_info / batch / "sample_info.csv"
    all_tasks = {
        "barcode_count": {},
        "escape_calc": {},
        "ref_merge": {},
        "sort_seq": {},
        "tite_seq": {},
    }
    all_df = []
    _df = (
        pd.read_csv(csvfile)
        .query("library in @libraries")
        .assign(batch=batch)
        .fillna("")
    )

    logger.info(f"Processing {batch}... Number of samples: {len(_df)}")
    if len(_df) == 0:
        logger.warning(f'No valid samples in batch {batch}')

    if len(_df) == 0:
        sys.stderr.write(f"Warning: {batch} - {csvfile}. No samples found.")
        return None, None

    if not necessary_columns.issubset(_df.columns):
        logger.error(
            f"Error: {batch} - {csvfile}. Missing necessary columns: {','.join(necessary_columns)}."
        )
        raise RuntimeError(f"Missing necessary columns in sample info.")

    if "fastq_files" not in _df.columns:  # try to find fq files
        is_cleaned = (raw / batch / "merge_fastq.py").exists()
        all_fq_files = []
        #all_candidate_files = [str(x) for x in (raw / batch).rglob(suffix if suffix[0] == "*" else "*" + suffix) if not '.low.f' in str(x)]
        all_candidate_files = [str(x) for x in glob.glob(str(raw/batch/'**'/(suffix if suffix[0] == "*" else "*" + suffix)), recursive=True) if not '.low.f' in str(x)]
        for sample in _df["sample"]:
            for k, v in sample_name_replace.items():
                sample = sample.replace(k, v)
            fq_files = get_fq_files(all_candidate_files, batch, sample, suffix, R2_pattern, is_cleaned)
            all_fq_files.append(fq_files)
        _df["fastq_files"] = all_fq_files
    
    # for col in ['description', 'AbConc']: # optional columns
    #     if col not in _df.columns:
    #         _df[col] = ''
    
    _df = _df.query("fastq_files != ''")  # remove samples without fastq files
    all_df.append(_df)

    if "group" not in _df.columns:
        _df = _df.assign(group=lambda x: x["library"])

    for (_lib, _grp), _samples in _df.query("is_ref == 'Y'").groupby(
        ["library", "group"]
    )["sample"]:
        all_tasks["ref_merge"][f"{batch}_{_lib}_{_grp}_REF_Merge"] = [
            batch + "_" + x for x in _samples
        ]

    for i in _df.index:
        row = _df.loc[i]
        if (row["is_ref"] == "N") and (
            not row["antibody"] in exclude_antibody_names
        ):
            _ref_name = f"{batch}_{row['library']}_{row['group']}_REF_Merge"
            if _ref_name in all_tasks["ref_merge"]:
                if batch not in all_tasks["escape_calc"]:
                    all_tasks["escape_calc"][batch] = {}
                all_tasks["escape_calc"][batch][batch + "_" + row["sample"]] = {
                    "library": row["library"],
                    "antibody": row["antibody"],
                    "ref": _ref_name,
                }
            else:
                logger.warning(
                    f"Warning: {batch} - {row['sample']}. Reference not found. Skipped."
                )

        all_tasks["barcode_count"][batch + "_" + row["sample"]] = {
            "library": row["library"],
            "batch": batch,
            "fastq_files": row["fastq_files"].split(","),
        }

    if (sample_info / batch / "sort_seq_info.csv").exists():
        all_tasks["sort_seq"].setdefault(batch, {})
        all_tasks["sort_seq"][batch].update(
            load_expression_info(
                sample_info / batch / "sort_seq_info.csv", prefix=f"{batch}_"
            )
        )

    if (sample_info / batch / "tite_seq_info.csv").exists():
        _libs = pd.unique(
            pd.read_csv(sample_info / batch / "tite_seq_info.csv")["library"]
        )
        for _lib in _libs:
            all_tasks["tite_seq"].append({"library": _lib, "batch": batch})

    return all_df, all_tasks

def main(config, batches, num_proc=1):
    logger = logging.getLogger(__name__)
    raw = Path(config["raw_bc"])
    sample_info = Path(config["sample_info_bc"])
    output = Path(config["output"])
    exclude_antibody_names = config["exclude_antibody_names"]

    libraries = list(config["libinfo"].keys())

    suffix = config["fq_suffix"] if "fq_suffix" in config else "*q.gz"
    R2_pattern = re.compile(
        r"((_2)|(_R2))\.f(ast)?q\.gz"
        if "R2_pattern" not in config
        else config["R2_pattern"]
    )

    logger.info(f"Libraries: {','.join(libraries)}")
    logger.info(f"fastq suffix: {suffix}")

    all_df = []
    all_tasks = {
        "barcode_count": {},
        "escape_calc": {},
        "ref_merge": {},
        "sort_seq": {},
        "tite_seq": {},
    }

    necessary_columns = set(["sample", "library", "antibody", "is_ref"])

    sample_name_replace = (
        config["sample_name_replace"]
        if "sample_name_replace" in config
        else {" ": "-", "/": "-", ".": "-"}
    )

    if num_proc == 1:
        for batch in batches:
            _df, _tasks = parse_batch(
                sample_info,
                raw,
                batch,
                necessary_columns,
                sample_name_replace,
                exclude_antibody_names,
                libraries,
                suffix,
                R2_pattern,
            )
            if _df is not None:
                all_df.append(_df)
                for k, v in _tasks.items():
                    all_tasks[k].update(v)
    else:
        with ProcessPoolExecutor(max_workers=num_proc) as executor:
            futures = []
            for batch in batches:
                futures.append(
                    executor.submit(
                        parse_batch,
                        sample_info,
                        raw,
                        batch,
                        necessary_columns,
                        sample_name_replace,
                        exclude_antibody_names,
                        libraries,
                        suffix,
                        R2_pattern,
                    )
                )

            for future in futures:
                _df, _tasks = future.result()
                if _df is not None:
                    all_df += _df
                    for k, v in _tasks.items():
                        all_tasks[k].update(v)

    yaml.dump(all_tasks, open(output / "_tasks.yaml", "w"))

    if len(all_df) > 0:
        pd.concat(all_df).to_csv(output / "sample_info_all.csv", index=None)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", help="config file")
    parser.add_argument("--workdir", help="working directory", default=".")
    parser.add_argument("--num_proc", help="number of processes", default=-1, type=int)

    args = parser.parse_args()
    config = yaml.safe_load(open(args.config, "r"))

    os.chdir(args.workdir)

    log_file = Path(config["output"]) / "logs" / "sample_info_clean.log.txt"
    log_file.parent.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(level=logging.INFO, filename=log_file, filemode="w")

    if "batches" in config and config["batches"] != "all":
        BATCHES = config["batches"].split(",")
    else:
        BATCHES = os.listdir(config["sample_info_bc"])

    logger = logging.getLogger(__name__)

    logger.info(f"Start: {time.ctime()}")
    logger.info(f"Sample Info Directory: {config['sample_info_bc']}")
    logger.info(f"Raw Barcode Directory: {config['raw_bc']}")
    logger.info(f"Output Directory: {config['output']}")

    logger.info(f"Batches: {BATCHES}")

    num_proc = int(args.num_proc if args.num_proc > 0 else (os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count())))
    logger.info(f"Number of processes: {num_proc}")

    main(config=config, batches=BATCHES, num_proc=num_proc)

    logger.info(f"End: {time.ctime()}")
    logging.shutdown()
