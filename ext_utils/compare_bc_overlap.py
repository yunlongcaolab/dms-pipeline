import gzip
import pandas as pd
from pathlib import Path

def read_and_count(fq_files: list[Path], length: int = 26) -> dict[str, int]:
    bc_dict = {}
    for fq_file in fq_files:
        assert fq_file.suffix == '.gz', f'{fq_file} is not a gzipped file'
        with gzip.open(fq_file, 'rt') as f:
            # parse fastq file
            line = f.readline()
            line = f.readline()
            while line:
                bc = line.strip()[:length]
                if bc not in bc_dict:
                    bc_dict[bc] = 1
                else:
                    bc_dict[bc] += 1
                for _ in range(3):
                    line = f.readline()
                line = f.readline()
    return bc_dict

def compare(bc_dict1: dict[str, int], bc_dict2: dict[str, int]):
    df1 = pd.DataFrame(bc_dict1.items(), columns=['barcode', 'count1'])
    df2 = pd.DataFrame(bc_dict2.items(), columns=['barcode', 'count2'])

    df = pd.merge(df1, df2, on='barcode', how='outer').fillna(0)

    print(len(bc_dict1), len(bc_dict2))
    print(f"Union: {len(df)}")
    print(f"Overlap: {len(df.query('count1 > 0 and count2 > 0'))}")

    df.to_csv('bc_overlap.csv', index=False)

def main():
    # read two lists of fq files
    fq1 = [f for f in Path("/lustre/grp/cyllab/share/DMS/TREM1/NGS_raw/20240621/N2413791_80-1551512759_2024-06-20/240618-A00151B").glob('*lib2*_1.f*q.gz')]
    # fq1 = [f for f in Path("/lustre/grp/cyllab/share/DMS/TREM1/NGS_raw/20240621/N2413791_80-1551512759_2024-06-20/240618-A00151B").glob('*lib1*_1.f*q.gz')]
    fq2 = [Path("/lustre/grp/cyllab/share/DMS/TREM1/NGS_raw/20240726/N2417131_80-1582634903_2024-07-25/240723-A00248A/T1-lib2_L1_1.fq.gz")]

    print(fq1)
    print(fq2)

    compare(
        read_and_count(fq1),
        read_and_count(fq2)
    )

if __name__ == '__main__':
    main()