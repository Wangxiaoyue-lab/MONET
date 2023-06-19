from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.gridspec import GridSpec

def show_colormaps(cmap_name=None, n=5, plot=True):
    """获得colormap基础颜色"""
    # 获取所有可用的颜色映射名称
    cmap_names = plt.colormaps()
    # 如果没有指定颜色映射名称，则展示所有颜色映射
    if cmap_name is None:
        cmap_name = cmap_names
    # 如果指定了颜色映射名称，则检查名称是否有效
    else:
        if isinstance(cmap_name, str):
            cmap_name = [cmap_name]
        for name in cmap_name:
            if name not in cmap_names:
                raise ValueError(f"Invalid colormap name: {name}")
    # 创建一个字典来存储颜色值
    colors_dict = {}
    # 获取颜色值
    for name in cmap_name:
        # 获取一个颜色映射对象
        cmap = plt.get_cmap(name)
        # 生成一些数值
        x = np.linspace(0, 1, n)
        # 将数值映射到颜色
        colors = cmap(x)
        # 将颜色值添加到字典中
        colors_dict[name] = colors
    # 如果 plot 参数为 True，则绘制条状图
    if plot:
        fig, axs = plt.subplots(len(cmap_name), 1, figsize=(n, len(cmap_name) * 0.5))
        axs = np.atleast_1d(axs)
        for ax, name in zip(axs, cmap_name):
            colors = colors_dict[name]
            # 绘制矩形并填充颜色
            for i, color in enumerate(colors):
                ax.add_patch(plt.Rectangle((i * 1.1 + 0.5, 0), 1, 0.4, color=color))
                if i == 0:
                    ax.text(i * 1.1 + 0.25, 0.2, name, ha="right", va="center")
            ax.set_xlim(0, n * 1.1 + 0.5)
            ax.set_ylim(0, 0.5)
            ax.axis("off")
        plt.show()
    # 返回颜色值字典
    return colors_dict



def show_markers(adata, keys_select: dict, marker_genes: dict, cluster: str, filename: str, n_rows: int = 3,n_cols: int = 3):
    with PdfPages(filename) as pdf:
        for ct in keys_select:
            print(f"{ct.upper()}:")
            fig, axs = plt.subplots(n_rows + 3, n_cols, figsize=(4*n_cols, 4*(n_rows + 5)), gridspec_kw={'wspace': 0.2, 'hspace': 0.2})
            
            # Create new axis in existing figure
            gs = GridSpec(n_rows + 3, n_cols, height_ratios=[4] + [0] + [0] + [1] * n_rows)
            ax2 = fig.add_subplot(gs[0, :])
            sc.pl.umap(data_sc, color=[cluster], legend_loc="on data", size=40, legend_fontsize = 20, ax=ax2, show=False)

            for i, gene in enumerate(marker_genes[ct]):
                ax = axs[(i // n_cols) + 3, i % n_cols]
                sc.pl.umap(
                    data_sc,
                    color=gene,
                    vmin=0,
                    vmax="p99",
                    sort_order=False,
                    frameon=False,
                    cmap="Reds",
                    ax=ax,
                    show=False
                )
                ax.set_title(gene)
            for i in range(len(marker_genes[ct]), n_rows*n_cols):
                ax = axs[(i // n_cols) + 3, i % n_cols]
                ax.axis('off')
            fig.subplots_adjust(left=0.04)
            fig.suptitle(ct.upper(), fontsize=60)
            pdf.savefig(fig)
            plt.close(fig)
            print("\n\n\n")



# https://scanpy-tutorials.readthedocs.io/en/latest/paga-paul15.html#Optional:-Denoising-the-graph
def show_paga_path(adata, celltype: str, lineages: dict, pseudotime: str, gene_names: List, filename: str):
    _, axs = pl.subplots(ncols=len(lineages), figsize=(6, 2.5), gridspec_kw={'wspace': 0.05, 'left': 0.12})
    pl.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)
    for ilineage, (descr, lineage) in enumerate(lineages):
        _, data = sc.pl.paga_path(
            adata, lineage, gene_names,
            show_node_names=False,
            ax=axs[ilineage],
            ytick_fontsize=12,
            left_margin=0.15,
            n_avg=50,
            annotations=[pseudotime],
            show_yticks=True if ilineage==0 else False,
            show_colorbar=False,
            color_map='Greys',
            groups_key=celltype,
            color_maps_annotations={pseudotime: 'viridis'},
            title='{} lineage'.format(descr),
            return_data=True,
            show=False)
        #data.to_csv('./write/paga_path_{}.csv'.format(descr))
    pl.savefig(filename)
    pl.show()

