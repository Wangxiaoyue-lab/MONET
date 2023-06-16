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
