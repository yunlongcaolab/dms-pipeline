import pandas as pd
import numpy as np
import os
from scipy.stats import norm, CensoredData
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

def _fit_log_normal(
    data_bins: dict[str, int],
    bins: dict[str, dict[str, float]],
    count_coef: int | None = None,
) -> tuple[str, float, float]:
    """Fit normal distribution to the log-transformed data.

    Args:
        data_bins (dict[str, int]): A dictionary mapping bin names to the number of cells in the bin
        bins (dict[str, dict[str, float]]): A dictionary mapping bin names to the left and right boundary of the bin
        count_coef (int, optional): Coefficient to multiply the count. Defaults to None (auto-determined based on the total counts).
    Returns:
        tuple[float, float]: A tuple of (mean, std) of the normal distribution
    """
    # construct censored data
    censored_data = []

    if count_coef is None:
        _total_counts = sum(data_bins.values())
        count_coef = min(100, 10000 / _total_counts)

    for bin, count in data_bins.items():
        left, right = bins[bin]["left"], bins[bin]["right"]
        left = 0 if left <= 0 else np.log(left)
        right = np.log(right)

        censored_data += int(count * count_coef + 0.5) * [[left, right]]

    censored_data = CensoredData(interval=censored_data)

    # fit normal distribution
    mean, std = norm.fit(censored_data)

    return mean, std


def load_expression_info(
    file: Path | str, prefix: str = ""
) -> dict[str, dict[str, dict[str, str | float | int]]]:
    """Load expression information from a csv file.

    Args:
        file (Path | str): Path to the csv file. The file should contain the following columns:
            - sample: Sample name
            - library: Library name
            - bin: bin of the sorting
            - num_cells: count of the cells in the bin
            - min_F: minimum fluorescence value of the bin
            - max_F: maximum fluorescence value of the bin
        prefix (str, optional): Prefix to add to the sample name. Defaults to ''.
    Returns:
        dict[str, dict[str, dict]]: A dictionary mapping library name to a dictionary that maps bin name to a dict of (sample, left, right, cells)
    """

    df = pd.read_csv(file)

    required_cols = ["sample", "library", "bin", "num_cells", "min_F", "max_F"]

    for col in required_cols:
        assert col in df.columns, f"Missing column {col}"

    lib_expr_info = {}
    for library, rows in df.groupby("library"):
        # bin id should be unique
        assert len(rows["bin"].unique()) == len(
            rows
        ), f"Duplicate bin id in library {library}"

        lib_expr_info[library] = {
            bin: {
                "sample": prefix + row["sample"],
                "min_F": float(row["min_F"]),
                "max_F": float(row["max_F"]),
                "n_cells": int(row["num_cells"]),
            }
            for bin, row in rows.set_index("bin").iterrows()
        }

    return lib_expr_info


class SortSeqData:
    """A class to store and analyze sorting data."""

    def __init__(
        self,
        bins: dict[str, float | int],
        bin_files: dict[str, str] = None,
        table: pd.DataFrame = None,
        merge_variant: bool = None,
    ):
        """Initialize the SortSeqData object."""
        self.bins = {}
        for bin, _info in bins.items():
            left = _info["min_F"]
            right = _info["max_F"]
            cells = _info["n_cells"]

            assert left < right, f"Invalid bin {bin}: left ({left}) >= right ({right})"
            assert left >= 0, f"Invalid bin {bin}: left = {left}"
            assert cells > 0, f"Invalid bin {bin}: cells = {cells}"

            self.bins[bin] = {
                "left": left,
                "right": right,
                "n_cells": cells,
            }

        if bin_files is not None:
            assert (
                merge_variant is not None
            ), "merge_variant should be specified when bin_files is provided"
            self.from_count_files(bin_files, table, merge_variant)

    def _normalize_count(self, df: pd.DataFrame) -> pd.DataFrame:
        """Normalize the count based on the total counts in each bin.

        Args:
            df (pd.DataFrame): A dataframe containing the following columns:
                - barcode: Barcode of the variant
                - aa_substitutions: Amino acid substitutions of the variant
                - bin: Bin
                - count: Count of reads of this barcode in the bin
        Returns:
            pd.DataFrame: A dataframe with the following columns added:
                - estimated_cells: Estimated number of cells in the bin
        """
        return (
            df.assign(
                bin_total=lambda x: x["bin"].map(
                    lambda b: self.bins[b]["total_counts"]
                ),
                bin_cells=lambda x: x["bin"].map(lambda b: self.bins[b]["n_cells"]),
            )
            .assign(
                estimated_cells=lambda x: x["count"] * x["bin_cells"] / x["bin_total"]
            )
            .drop(columns=["bin_total", "bin_cells"])
        )

    def from_count_files(
        self,
        files: dict[str, Path | str],
        table: pd.DataFrame = None,
        merge_variant: bool = True,
    ):
        """Load data from csv files.

        Args:
            files (dict[str, Path]): A dictionary mapping bin to the path of the csv file. Each csv file should contain the following two columns:
                - barcode: Barcode of the variant
                - count: Count of reads of this barcode in the bin
            table (pd.DataFrame, optional): A table containing the following columns:
                - barcode: Barcode of the variant
                - aa_substitutions: Amino acid substitutions of the variant
                If None, the files should contain the aa_substitutions column. Defaults to None.
            merge_variant (bool, optional): Whether to merge the same variants with different barcodes. Defaults to True
        """

        df = []

        for bin, file in files.items():
            df.append(
                pd.read_csv(file)
                .assign(bin=bin)
                .merge(table[["barcode", "aa_substitutions"]], on="barcode", how="left")
            )

        self.from_df(pd.concat(df), merge_variant)

    def from_df(self, df: pd.DataFrame, merge_variant: bool = True):
        """Load data from a dataframe.

        Args:
            df (pd.DataFrame): A dataframe containing the following columns:
                - barcode: Barcode of the variant
                - aa_substitutions: Amino acid substitutions of the variant
                - bin: Bin of the sorting
                - count: Count of reads of this barcode in the bin
            merge_variant (bool, optional): Whether to merge the same variants with different barcodes. Defaults to True
        """
        required_cols = ["barcode", "aa_substitutions", "bin", "count"]

        for col in required_cols:
            assert col in df.columns, f"Missing column {col}"

        df["aa_substitutions"] = df["aa_substitutions"].fillna("")
        df["count"] = df["count"].fillna(0)

        if merge_variant:
            df = (
                df.groupby(["aa_substitutions", "bin"])["count"]
                .sum()
                .reset_index()
                .assign(barcode=lambda x: x["aa_substitutions"].str.replace(" ", "+"))
            )

        # add bin total counts
        for bin, value in df.groupby("bin")["count"].sum().to_dict().items():
            self.bins[bin]["total_counts"] = value
        self.df = self._normalize_count(df).sort_values(
            "estimated_cells", ascending=False
        )

    def fit_log_normal(self, num_proc: int = 1):
        """Log-transform and fit normal distribution to all variants data saved in the object DataFrame.

        Args:
            num_proc (int, optional): Number of processes to use. Defaults to 1.
        """
        num_proc = min(int(num_proc), len(self.df["barcode"].unique()), os.cpu_count())
        if num_proc > 1:
            with ProcessPoolExecutor(num_proc) as p:
                futures = []
                for barcode, rows in self.df.groupby("barcode"):
                    data_bins = rows.set_index("bin")["estimated_cells"].to_dict()
                    futures.append(
                        (barcode, p.submit(_fit_log_normal, data_bins, self.bins))
                    )
                results = {barcode: future.result() for barcode, future in futures}
        else:
            results = {}
            for barcode, rows in self.df.groupby("barcode"):
                data_bins = rows.set_index("bin")["estimated_cells"].to_dict()
                results[barcode] = _fit_log_normal(data_bins, self.bins)

        self.df_fit = (
            self.df.pivot(
                index=["barcode", "aa_substitutions"],
                columns="bin",
                values="estimated_cells",
            )
            .reset_index()
            .merge(
                pd.DataFrame(results, index=["mean", "std"]).T.reset_index(
                    names="barcode"
                ),
                on="barcode",
                how="left",
            )
        )


def main(snakemake):
    input_files = snakemake.input

    print(input_files)
    print(snakemake.params)
    print(snakemake.params.bins)

    bin_files = {}
    for bin in snakemake.params.bins:
        for file in input_files:
            if bin in file:
                bin_files[bin] = file
                break
        else:
            raise ValueError(f"Bin {bin} not found in input")

    sort_seq_data = SortSeqData(
        bins=snakemake.params.bins,
        bin_files=bin_files,
        table=pd.read_csv(snakemake.params.table),
        merge_variant=snakemake.config["merge_variant"],
    )

    sort_seq_data.fit_log_normal(num_proc=snakemake.resources.cpu_per_task)

    if snakemake.config["merge_variant"]:
        # sort_seq_data.df.drop(columns='barcode').pivot(index='aa_substitutions', columns='bin', values='estimated_cells').reset_index().to_csv(snakemake.output.counts, index=False)
        sort_seq_data.df_fit.drop(columns="barcode").to_csv(
            snakemake.output.variant, index=False, float_format="%.7g"
        )
    else:
        # sort_seq_data.df.pivot(index=['barcode','aa_substitutions'], columns='bin', values='estimated_cells').reset_index().to_csv(snakemake.output.counts, index=False)
        sort_seq_data.df_fit.to_csv(
            snakemake.output.variant, index=False, float_format="%.7g"
        )


if __name__ == "__main__":
    main(snakemake)
