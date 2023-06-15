from scipy.stats import median_abs_deviation


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
    print(f"Lower bound: {lower_bound}")
    print(f"Upper bound: {upper_bound}")
    return outlier
