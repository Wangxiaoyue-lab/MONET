from typing import Optional
import numpy as np
from scipy.stats import median_abs_deviation


def check_outlier(adata, metric: Optional[str] = None, nmads: int = 5, mt_lower: int = 5, bound: str = "both"):
    def is_outlier(adata, metric: str, nmads: int, bound: str = "both"):
        M = adata.obs[metric]
        lower_bound = np.median(M) - nmads * median_abs_deviation(M)
        upper_bound = np.median(M) + nmads * median_abs_deviation(M)
        if bound == "lower":
            outlier = M < lower_bound
        elif bound == "upper":
            outlier = upper_bound < M
        else:
            outlier = (M < lower_bound) | (upper_bound < M)
        print(f"{metric}--Lower bound: {lower_bound}")
        print(f"{metric}--Upper bound: {upper_bound}")
        return outlier

    if metric is None:
        adata.obs["counts_coutlier"] = (
            is_outlier(adata, "log1p_total_counts", nmads, bound = "upper")
            | is_outlier(adata, "log1p_n_genes_by_counts", nmads, bound = "upper")
            | is_outlier(adata, "pct_counts_in_top_20_genes", nmads, bound = "upper")
        )
        adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3, bound="upper") | (
            adata.obs["pct_counts_mt"] > mt_lower
        )
    else:
        adata.obs[f"{metric}_outlier"] = is_outlier(adata, metric, nmads, bound: str = "both")

    return adata


def check_marker(adata, markers_genes: dict):
    marker_genes_in_data = dict()
    for ct, markers in marker_genes.items():
        markers_found = list()
        for marker in markers:
            if marker in adata.var.index:
                markers_found.append(marker)
        marker_genes_in_data[ct] = markers_found

    return marker_genes_in_data