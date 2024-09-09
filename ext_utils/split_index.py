import argparse
import gzip
from pathlib import Path

parser = argparse.ArgumentParser(description='Split unsplitted fastq file into files according to an index dictionary')

parser.add_argument('fastq', type=str, help='Input fastq file')
parser.add_argument('index', type=str, help='Index dictionary with two columns split by whitespace (sample, index)')
parser.add_argument('output', type=str, help='Output directory')
parser.add_argument('--index-length', type=int, default=None, help='Length of index sequence')
parser.add_argument('--index-pos', type=str, default="end", help='Position of index sequence: start or end')
parser.add_argument('--no-rev-comp', action='store_true', help='Do not reverse complement index sequence')


def reverse_complement(seq):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement[base] for base in reversed(seq))


def main(args):
    indices = {}
    blacklist_index = set()
    blacklist_samples = set()
    with open(args.index, "r", encoding="utf-8") as f:
        for line in f:
            sample, index = line.strip().split()
            index = (
                index.strip() if args.no_rev_comp else reverse_complement(index.strip())
            )
            index = index[0 : args.index_length] if args.index_length else index

            if index in indices:
                blacklist_index.add(index)
                blacklist_samples.add(indices[index])
                blacklist_samples.add(sample)
            else:
                indices[index] = sample

    if len(blacklist_index) > 0 or len(blacklist_samples) > 0:
        print("Ambiguous indices:")
        for index in blacklist_index:
            print(f"  {index}")
            indices.pop(index)
        print("Invalid samples:")
        for sample in blacklist_samples:
            print(f"  {sample}")

    # all indices of the same length
    assert len(set(indices.values())) == len(
        indices
    ), "All samples should have unique indices"

    # all sample names should be unique
    assert len(set(indices.values())) == len(
        set(indices.keys())
    ), "All samples should have unique names"

    print("Using index dictionary:")
    for index, sample in indices.items():
        print(f"  {index} -> {sample}")

    args.index_length = (
        len(next(iter(indices.keys()))) if not args.index_length else args.index_length
    )

    # open output files
    out_files = {}
    Path(args.output).mkdir(parents=True, exist_ok=True)
    for sample in set(indices.values()):
        out_files[sample] = gzip.open(
            Path(args.output) / f"{sample}.fastq.gz", "wt", encoding="utf-8"
        )

    out_files["unknown"] = gzip.open(
        Path(args.output) / "unknown.fastq.gz", "wt", encoding="utf-8"
    )

    if Path(args.fastq).suffix == ".gz":
        open_func = gzip.open
    else:
        open_func = open

    with open_func(args.fastq, "rt", encoding="utf-8") as f:
        for line in f:
            if line.startswith("@"):
                header = line
                seq = next(f).strip()
                plus = next(f)
                qual = next(f)
                _index = (
                    seq[-args.index_length :]
                    if args.index_pos == "end"
                    else seq[: args.index_length]
                )
                sample = indices.get(_index, "unknown")

                if _index in indices:
                    seq = (
                        seq[: -args.index_length]
                        if args.index_pos == "end"
                        else seq[args.index_length :]
                    )

                out_files[sample].write(f"{header}{seq}\n{plus}{qual}")

    for f in out_files.values():
        f.close()


if __name__ == "__main__":
    main(parser.parse_args())
