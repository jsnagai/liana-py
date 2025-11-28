import anndata
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy 
import itertools
from liana._constants import DefaultValues as V
from liana._constants import Keys as K
from liana._docs import d

from liana._logging import _logg

@d.dedent
def spatial_inflow(
    adata: anndata.AnnData = None,
    groupby: str = None,
    spatial_key = K.spatial_key,
    labels: list[str] = None, 
    interactions = K.interactions,
    figure_size: tuple = (10, 6),
    normalize: bool = True,
    percentile_scaling: tuple[int, int] | None = None, 
    show_counts: bool = True,
):

    """
    Plot inflow scores for single interaction across spatial coordinates.

    Parameters
    ---------- 
    %(adata)s
    %(groupby)s
    spot_size
        Size of the spots in the scatter plot
    labels
        List of labels to compare
    %(spatial_key)s
    %(interactions)s
    %(cmap)s
    %(figure_size)s
    normalize
        Whether to normalize expression values for color mapping
    percentile_scaling
        Tuple specifying percentiles for scaling color mapping
    show_counts
        Whether to display counts of cells per label.
    """

    # Validate inputs
    if not labels:
        raise ValueError("labels must contain at least one label")
    if spatial_key not in adata.obsm:
        raise KeyError(f"'{spatial_key}' not found in adata.obsm")

    # Default colormaps if not provided
    default_cmaps = ['Blues', 'Reds', 'Greens', 'Purples', 'Oranges',
                            'YlOrBr', 'PuRd', 'BuGn', 'GnBu', 'OrRd']
    cmaps = list(itertools.islice(itertools.cycle(default_cmaps), len(labels)))


    # Prepare data
    coords = adata.obsm[spatial_key]
    cell_type_data = []
    
    for label in labels:
        mask = adata.obs[groupby].eq(label)
        if not mask.any():
            _logg(f"No valid label found for {label}", stacklevel=2)
            continue

        coords_type = coords[mask.values, :]

        if interactions in adata.var_names:
            expr_type = adata[mask][:, interactions].X.toarray().ravel()
        else:
            _logg(f"'{interactions}' not found in adata.var_names. Skipping {label}.", stacklevel=2)
            continue

        cell_type_data.append({
            'label': label,
            'coords': coords_type,
            'expression': expr_type.copy(),
            'count': (expr_type > 0).sum()
        })

    # Normalize and scale expression data
    for data in cell_type_data:
        expr = data['expression']
        # Apply percentile clipping
        if percentile_scaling is not None:
            low_p, high_p = percentile_scaling
            scale_min = numpy.percentile(expr, low_p)
            scale_max = numpy.percentile(expr, high_p)
            expr = numpy.clip(expr, scale_min, scale_max)
        else:
            scale_min = numpy.min(expr)
            scale_max = numpy.max(expr)

        # Normalize to 0-1 
        if normalize:
            if scale_max - scale_min > 0:
                expr = (expr - scale_min) / (scale_max - scale_min)
            else:
                expr = numpy.zeros_like(expr)

        data['expression'] = expr

    fig = plt.figure(figsize=figure_size)
    ax = fig.add_axes([0.1, 0.24, 0.8, 0.72])
    ax.scatter(coords[:, 0], coords[:, 1],
               color='lightgrey', s=3, alpha=0.2, rasterized=True)
    # 2. Scatter plots for each cell type
    scatter_objects = []
    for i, data in enumerate(cell_type_data):
        sc = ax.scatter(
            data['coords'][:, 0], data['coords'][:, 1],
            c=data['expression'], cmap=cmaps[i], s=3,
            label=data['label'], alpha=0.8,
            rasterized=True,
        )
        scatter_objects.append(sc)

    n_bars = len(scatter_objects)
    total_width = 0.8
    bar_margin = 0.02 # space between bars
    bar_width = (total_width - (n_bars - 1) * bar_margin) / n_bars
    for i, (sc, data) in enumerate(zip(scatter_objects, cell_type_data, strict=True)):
        left = 0.1 + i * (bar_width + bar_margin)
        cax = fig.add_axes([left, 0.08, bar_width, 0.03])
        cb = plt.colorbar(sc, cax=cax, orientation='horizontal')

        label_text = data['label']
        if show_counts:
            label_text += f" (n={data['count']})"

        cb.set_label(label_text, fontsize=9)

        cb.ax.tick_params(labelsize=8)

    title_text = f"{interactions} inflow"
    ax.set_title(title_text, fontsize=14, pad=10)
    
    return fig, ax