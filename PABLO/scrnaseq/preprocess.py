from typing import Optional
import anndata2ri
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
import numpy as np
from scipy.stats import median_abs_deviation

def mark_outlier(adata, metric: Optional[str] = None, nmads: int = 5, mt_lower: int = 5, bound: str = "both"):
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

def run_norm2scale(adata, layer_raw: str = None,n_top_genes: int = 2000, regress_vars: List = None, **kwargs):
    if layer_raw:
        adata.layers['raw_counts'] = adata.X
    else:
        adata.X = adata.layers[layer_raw]

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["log1p_norm"] = adata.X
    sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes = n_top_genes)
    adata = adata[:, adata.var.highly_variable]
    if regress_vars:
        regress_vars = ['total_counts', 'pct_counts_mt']
    sc.pp.regress_out(adata, regress_vars)
    sc.pp.scale(adata, max_value=10)
    adata.layers["scale"] = adata.X
    return adata

def run_sctransform(adata, layer=None, **kwargs):
    """
    Run SCTransform on the adata object and add the results as a new layer
    Arguments:
    adata: AnnData object containing UMI counts matrix
    layer: Name of the layer to which the results will be added; default is None
    **kwargs: Additional parameters passed to sctransform

    Keyword arguments:
    cell.attr: A metadata with cell attributes
    reference.SCT.model: If not NULL, compute residuals for the object using the provided SCT model; supports only log_umi as the latent variable. If residual.features are not specified, compute for the top variable.features.n specified in the model which are also present in the object. If residual.features are specified, the variable features of the resulting SCT assay are set to the top variable.features.n in the model.
    do.correct.umi: Place corrected UMI matrix in assay counts slot; default is TRUE
    ncells: Number of subsampling cells used to build NB regression; default is 5000
    residual.features: Genes to calculate residual features for; default is NULL (all genes). If specified, will be set to VariableFeatures of the returned object.
    variable.features.n: Use this many features as variable features after ranking by residual variance; default is 3000. Only applied if residual.features is not set.
    variable.features.rv.th: Instead of setting a fixed number of variable features, use this residual variance cutoff; this is only used when variable.features.n is set to NULL; default is 1.3. Only applied if residual.features is not set.
    vars.to.regress: Variables to regress out in a second non-regularized linear regression. For example, percent.mito. Default is NULL
    do.scale: Whether to scale residuals to have unit variance; default is FALSE
    do.center: Whether to center residuals to have mean zero; default is TRUE
    clip.range: Range to clip the residuals to; default is c(-sqrt(n/30), sqrt(n/30)), where n is the number of cells
    vst.flavor: When set to 'v2' sets method = glmGamPoi_offset, n_cells=2000, and exclude_poisson = TRUE which causes the model to learn theta and intercept only besides excluding poisson genes from learning and regularization
    conserve.memory: If set to TRUE the residual matrix for all genes is never created in full; useful for large data sets, but will take longer to run; this will also set return.only.var.genes to TRUE; default is FALSE
    return.only.var.genes: If set to TRUE the scale.data matrices in output assay are subset to contain only the variable genes; default is TRUE
    seed.use: Set a random seed. By default, sets the seed to 1448145. Setting NULL will not set a seed.
    verbose: Whether to print messages and progress bars
    assay: Name of assay to pull the count data from; default is 'RNA'
    new.assay.name: Name for the new assay containing the normalized data; default is 'SCT'
    """
    # Extract the representation matrix layer from the anndata object
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
    r("seurat_obj@assays$SCT@var.features")

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

    # Add the SCT variance to the anndata object
    adata.var["sct_hvg"] = np.isin(
        np.array(adata.var_names), np.array(r("seurat_obj@assays$SCT@var.features"))
    )
    return adata
