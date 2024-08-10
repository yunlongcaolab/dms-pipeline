import pandas as pd
import numpy as np
import scipy.stats

def _fit_log_normal(barcode, data_bins: dict[str, int], bins: dict[str, dict[str, float]]):
    # construct censored data
    censored_data = []
    
    for bin, count in data_bins.items():
        left, right = bins[bin]['left'], bins[bin]['right']
        left = -np.inf if left == 0 else np.log10(left)
        right = np.log10(right)

        censored_data += count * [[left, right]]

    censored_data = scipy.stats.CensoredData(interval=censored_data)

    # fit normal distribution
    mean, std = scipy.stats.norm.fit(censored_data)

    return barcode, mean, std

class SortSeqData:
    def __init__(self, bins: dict[str, tuple[int, int]]):
        self.bins = {}
        for bin, (left, right, cells) in bins.items():
            assert left < right, f"Invalid bin {bin}: left ({left}) >= right ({right})"
            assert left >= 0, f"Invalid bin {bin}: left = {left}"
            assert cells > 0, f"Invalid bin {bin}: cells = {cells}"

            self.bins[bin] = {
                'left': left,
                'right': right,
                'n_cells': cells,
            }
    
    def _normalize_count(self, df: pd.DataFrame):
        return df.assign(
            bin_total=lambda x: x['bin'].map(lambda b: self.bins[b]['total_counts']),
            bin_cells=lambda x: x['bin'].map(lambda b: self.bins[b]['n_cells']),
        ).assign(
            estimated_cells=lambda x: (x['count'] * x['bin_cells'] / x['bin_total'] + 0.5).astype(int)
        ).drop(columns=['bin_total', 'bin_cells'])

    def from_df(self, df: pd.DataFrame, merge_variant: bool = True):
        required_cols = ['barcode', 'aa_substitutions', 'bin', 'count']

        for col in required_cols:
            assert col in df.columns, f"Missing column {col}"
        
        df['aa_substitutions'] = df['aa_substitutions'].fillna('')
        df['count'] = df['count'].fillna(0)

        if merge_variant:
            df = df.groupby(['aa_substitutions', 'bin']).agg({
                'barcode': lambda x: ','.join(x),
                'count': 'sum'
            }).reset_index()
        
        # add bin total counts
        for bin, value in df.groupby('bin')['count'].sum().to_dict().items():
            self.bins[bin]['total_counts'] = value
        self.df = self._normalize_count(df)

    def fit_log_normal(self, num_proc: int = 1):
        if num_proc > 1:
            from multiprocessing import Pool
            params = []
            with Pool(num_proc) as p:
                for barcode, rows in self.df.groupby('barcode'):
                    data_bins = rows.set_index('bin')['estimated_cells'].to_dict()
                    params.append((barcode, data_bins, self.bins))
                chunksize = len(params) // num_proc + 1
                results = p.starmap(_fit_log_normal, params, chunksize=chunksize)
        else:
            results = []
            for barcode, rows in self.df.groupby('barcode'):
                data_bins = rows.set_index('bin')['estimated_cells'].to_dict()
                results.append(_fit_log_normal(barcode, data_bins, self.bins))
        
        self.df_fit = self.df.pivot(index=['barcode', 'aa_substitutions'], columns='bin', values='estimated_cells').reset_index().merge(
            pd.DataFrame(results, columns=['barcode', 'mean', 'std']),
            on='barcode',
            how='left'
        )

def main():

    pass

if __name__ == "__main__":
    main()