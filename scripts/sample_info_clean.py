import pandas as pd
import sys, os, re
import yaml
import argparse
import time

from pathlib import Path
from sort_seq import load_expression_info

parser = argparse.ArgumentParser()
parser.add_argument("--config", help="config file")
parser.add_argument("--workdir", help="working directory", default=".")


def main(config, batches, logobj):
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

    logobj.write(f"Libraries: {','.join(libraries)}\n")
    logobj.write(f"fastq suffix: {suffix}\n")
    logobj.flush()

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

    for batch in batches:
        csvfile = sample_info / batch / "sample_info.csv"
        _df = (
            pd.read_csv(csvfile)
            .query("library in @libraries")
            .assign(batch=batch)
            .fillna("")
        )

        print(f"Processing {batch}... Number of samples: {len(_df)}")

        if len(_df) == 0:
            sys.stderr.write(f"Warning: {batch} - {csvfile}. No samples found.\n")
            continue

        if not necessary_columns.issubset(_df.columns):
            logobj.write(
                f"Error: {batch} - {csvfile}. Missing necessary columns: {','.join(necessary_columns)}.\n"
            )
            logobj.flush()
            raise RuntimeError(f"Missing necessary columns in sample info.")

        if "fastq_files" not in _df.columns:  # try to find fq files
            all_fq_files = []
            for sample in _df["sample"]:
                for k, v in sample_name_replace.items():
                    sample = sample.replace(k, v)
                fq_files = [str(x) for x in (raw / batch).rglob(f"*{sample}" + suffix)]
                if len(fq_files) == 0:
                    logobj.write(
                        f"Warning: {batch} - {sample}. FASTQ not found. Try alternative path. \n"
                    )
                    fq_files = [
                        str(x)
                        for x in (raw / batch).rglob(
                            os.path.join(
                                f"*{sample}*",
                                suffix if suffix[0] == "*" else "*" + suffix,
                            )
                        )
                    ]

                if len(fq_files) == 0:
                    logobj.write(
                        f"Warning: {batch} - {sample}. FASTQ not found in alternative path. Skipped. \n"
                    )
                    all_fq_files.append("")
                elif "," in "".join(fq_files):
                    logobj.write(
                        f"Error: {batch} - {sample}. Comma in fastq file path.\n"
                    )
                    raise RuntimeError(f"Comma in fastq file path: {batch} - {sample}")
                elif len(fq_files) > 1:
                    logobj.write(
                        f"Warning: {batch} - {sample}. Multiple FASTQ found. {','.join(fq_files)}. Try removing R2 files.\n"
                    )
                    fq_files = [x for x in fq_files if not R2_pattern.search(x)]

                    logobj.write(
                        f"    {len(fq_files)} files found after removing. {','.join(fq_files)}.\n"
                    )
                    all_fq_files.append(",".join(fq_files))
                else:
                    all_fq_files.append(fq_files[0])
                logobj.flush()
            _df["fastq_files"] = all_fq_files
        # for col in ['description', 'AbConc']: # optional columns
        #     if col not in _df.columns:
        #         _df[col] = ''
        _df = _df.query("fastq_files != ''")  # remove samples without fastq files
        all_df.append(
            # _df[['batch', 'sample', 'library', 'antibody', 'is_ref', 'description', 'AbConc', 'fastq_files']]
            _df
        )

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
                    logobj.write(
                        f"Warning: {batch} - {row['sample']}. Reference not found. Skipped.\n"
                    )
                    logobj.flush()

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
    yaml.dump(all_tasks, open(output / "_tasks.yaml", "w"))

    if len(all_df) > 0:
        pd.concat(all_df).to_csv(output / "sample_info_all.csv", index=None)


if __name__ == "__main__":
    args = parser.parse_args()
    config = yaml.safe_load(open(args.config, "r"))

    os.chdir(args.workdir)

    log_file = Path(config["output"]) / "logs" / "sample_info_clean.log.txt"
    log_file.parent.mkdir(parents=True, exist_ok=True)

    if "batches" in config and config["batches"] != "all":
        BATCHES = config["batches"].split(",")
    else:
        BATCHES = os.listdir(config["sample_info_bc"])

    with open(log_file, "w") as logobj:

        logobj.write(f"Start: {time.ctime()}\n")
        logobj.write(f"Sample Info Directory: {config['sample_info_bc']}\n")
        logobj.write(f"Raw Barcode Directory: {config['raw_bc']}\n")
        logobj.write(f"Output Directory: {config['output']}\n")

        logobj.write(f"Batches: {BATCHES}\n")

        main(config=config, batches=BATCHES, logobj=logobj)
        logobj.write(f"End: {time.ctime()}\n")
