#def ad2seurat:


#def seurat2ad:

def obj_ad2sce(adata):
    """transfer anndata to sce(R)"""
    import anndata2ri
    from rpy2.robjects import r
    anndata2ri.activate()
    r.assign('adata',adata)
    sce = r('as(adata, "SingleCellExperiment")')
    return sce

def obj_sce2ad(sce):
    """transfer sce(R) to anndata"""
    import anndata2ri
    from rpy2.robjects import r
    anndata2ri.activate()
    with anndata2ri.convertrt:
        adata = ro.conversion.rpy2py(sce)
    return adata


def obj_sce2seurat(sce):
    from rpy2.robjects import r
    from rpy2.robjects.packages import importr
    seurat = importr('Seurat')
    obj_seurat = seurat.as_Seurat(sce)
    return obj_seurat