import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np


def split_umap(adata, split_by, ncol=2, nrow=None, **kwargs):
    categories = adata.obs[split_by].cat.categories
    if nrow is None:
        nrow = int(np.ceil(len(categories) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(5 * ncol, 4 * nrow))
    axs = axs.flatten()
    for i, cat in enumerate(categories):
        ax = axs[i]
        sc.pl.umap(
            adata[adata.obs[split_by] == cat], ax=ax, show=False, title=cat, **kwargs
        )
    plt.tight_layout()
