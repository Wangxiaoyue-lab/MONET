
def mark_outlier(adata, metric: str, nmads: int):
    from scipy.stats import median_abs_deviation
    def is_outlier(adata, metric: str, nmads: int):
        M = adata.obs[metric]
        outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
            np.median(M) + nmads * median_abs_deviation(M) < M
        )
        return outlier

    adata.obs["counts_coutlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
        adata.obs["pct_counts_mt"] > 8
    )
    return adata



def run_sctransform(adata, layer=None, **kwargs):
    # Extract the representation matrix layer from the anndata object
    import anndata2ri
    from rpy2.robjects.packages import importr
    from rpy2.robjects import r, pandas2ri
    import numpy as np

    anndata2ri.activate()
    pandas2ri.activate()
    if layer:
        mat = adata.layers[layer]
    else:
        mat = adata.X

    # Set the column names of the input matrix
    cell_names = adata.obs_names
    gene_names = adata.var_names
    r.assign("mat", mat.T)
    r.assign("cell_names", cell_names)
    r.assign("gene_names", gene_names)
    r("colnames(mat) <- cell_names")
    r("rownames(mat) <- gene_names")

    # Convert the matrix into a Seurat object
    seurat = importr("Seurat")
    r("seurat_obj <- CreateSeuratObject(mat)")

    # Run SCTransform
    for k, v in kwargs.items():
        r.assign(k, v)
    kwargs_str = ", ".join([f"{k}={k}" for k in kwargs.keys()])
    r(f'seurat_obj <- SCTransform(seurat_obj,vst.flavor="v2", {kwargs_str})')

    # Extract the SCT data and add it as a new layer in the original anndata object
    sct_data = np.asarray(r["as.matrix"](r("seurat_obj@assays$SCT@data")))
    adata.layers["SCT_data"] = sct_data.T
    sct_data = np.asarray(r["as.matrix"](r("seurat_obj@assays$SCT@counts")))
    adata.layers["SCT_counts"] = sct_data.T
    return adata
