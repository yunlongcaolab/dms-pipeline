import pandas as pd
import glob, os, sys
import yaml
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--config", help="config file")

def main(raw, sample_info, batches, output, logobj, exclude_antibody_names=[]):
    all_df = []
    all_tasks = {
        "barcode_count": {},
        "escape_calc": {},
        "ref_merge": {},
        "sort_seq": [],
        "tite_seq": [],
    }

    for batch in batches:
        csvfile = os.path.join(sample_info, batch, 'sample_info.csv')
        _df = pd.read_csv(csvfile).assign(batch=batch)
        if 'fastq_files' not in _df.columns: # try to find fq files
            all_fq_files = []
            for sample in _df['sample']:
                fq_files = glob.glob(os.path.join(raw, '**',f'*{sample}*q.gz'), recursive=True)
                if (len(fq_files) == 0):
                    logobj.write(f"Warning: {batch} - {sample}. FASTQ not found. Skipped. \n")
                    all_fq_files.append('')
                elif ',' in ''.join(fq_files):
                    logobj.write(f"Error: {batch} - {sample}. Comma in fastq file path.\n")
                    raise RuntimeError(f"Comma in fastq file path: {batch} - {sample}")
                elif (len(fq_files) > 1):
                    logobj.write(f"Warning: {batch} - {sample}. Multiple FASTQ found. {','.join(fq_files)}\n")
                    all_fq_files.append(','.join(fq_files))
                else:
                    all_fq_files.append(fq_files[0])
            _df['fastq_files'] = all_fq_files

        for col in ['description', 'AbConc']: # optional columns
            if col not in _df.columns:
                _df[col] = ''
        _df = _df.query("fastq_files != ''") # remove samples without fastq files
        all_df.append(
            _df[['batch', 'sample', 'library', 'antibody', 'is_ref', 'description', 'AbConc', 'fastq_files']]
        )

        for i in _df.index:
            row = _df.loc[i]
            if (row['is_ref'] == "N") and (not row['antibody'] in exclude_antibody_names):
                if batch not in all_tasks['escape_calc']:
                    all_tasks['escape_calc'][batch] = {}
                all_tasks['escape_calc'][batch][batch+'_'+row["sample"]] = {
                        "library": row['library'],
                        "antibody": row['antibody'],
                        "ref": batch+f"_{row['library']}_REF_Merge"
                    }
            all_tasks['barcode_count'][batch+"_"+row["sample"]] = {
                    "library": row['library'],
                    "batch": batch,
                    "fastq_files": row['fastq_files'].split(','),
                }

        for _refs in _df.query("is_ref == 'Y'").groupby('library')['sample']:
            all_tasks['ref_merge'][batch+"_"+_refs[0]+"_REF_Merge"] = [batch+'_'+x for x in _refs[1]]

        if os.path.exists(os.path.join(sample_info, batch, 'sort_seq.csv')):
            _libs = pd.unique(pd.read_csv(os.path.join(sample_info, batch, 'sort_seq_info.csv'))['library'])
            for _lib in _libs:
                all_tasks['sort_seq'].append(
                    {
                        "library": _lib,
                        "batch": batch
                    }
                )
        
        if os.path.exists(os.path.join(sample_info, batch, 'tite_seq.csv')):
            _libs = pd.unique(pd.read_csv(os.path.join(sample_info, batch, 'tite_seq_info.csv'))['library'])
            for _lib in _libs:
                all_tasks['tite_seq'].append(
                    {
                        "library": _lib,
                        "batch": batch
                    }
                )
    yaml.dump(all_tasks, open(os.path.join(output, "_tasks.yaml"), 'w'))
    pd.concat(all_df).to_csv(os.path.join(output, 'sample_info_all.csv'), index=None)

if __name__ == "__main__":
    args = parser.parse_args()
    config = yaml.safe_load(open(args.config, 'r'))
    BATCHES = os.listdir(config['sample_info_bc'])

    log_file = os.path.join(config['output'], "logs", "sample_info_clean.log.txt")
    if not os.path.exists(os.path.dirname(log_file)):
        os.makedirs(os.path.dirname(log_file))

    with open(log_file, 'w') as logobj:
        main(
            raw = config['raw_bc'], 
            sample_info = config['sample_info_bc'],
            batches = BATCHES,
            output=config['output'], 
            logobj=logobj,
            exclude_antibody_names=config['exclude_antibody_names']
        )