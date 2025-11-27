from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from anndata import AnnData 

from liana._constants import DefaultValues as V
from liana._constants import Keys as K
from liana._docs import d


@d.dedent
def lr_heatmap(
    adata: AnnData = None,
    ligand: str = None,
    receptor: str = None,
    groupby:  str = None,
    lr_sep: str = V.lr_sep,
    clip: bool = True,
    filter_min_mean: float | None = None,
    title: str | None = None,
    xlabel: str = "Receiver cell type",
    ylabel: str = "Sender cell type",
    cmap: str = V.cmap,
    figure_size: tuple = (5, 5),
    filter_fun: callable = None,
    **kwargs: Any
):
    """
    Heatmap of ligand-receptor (LR) interaction by source and target cell types

    Parameters
    ----------
    %(adata)s
    ligand : str
        Name of the ligand gene.
    receptor : str
        Name of the receptor gene.
    %(groupby)s
    lr_sep: str
            Separator to use for interaction names.
    clip : bool
        Whether to clip the color scale to the 5th and 95th percentiles.
    filter_min_mean : float | None, optional
        Minimum mean value across rows/columns to retain them in the heatmap. If None, no filtering is applied. Defaults to None.
    
    title : str | None, optional
        Title of the heatmap. If None, a default title is used. Defaults to None.
    xlabel : str, optional
        Label for the x-axis. Defaults to "Receiver cell type".
    ylabel : str, optional
        Label for the y-axis. Defaults to "Sender cell type".

    %(filter_fun)s
    %(cmap)s
    %(figure_size)s

    Returns
    -------
        pd.DataFrame
            The processed heatmap data used for plotting.

        A `plotnine.ggplot` instance
    """

    # 1. Input Validation and Data Subsetting
    
    if not isinstance(adata, AnnData):
        raise TypeError("adata must be an AnnData object.")
    
    if groupby not in adata.obs.columns:
        raise ValueError(f"'{groupby}' not found in adata.obs. Available columns are: {list(adata.obs.columns)}")
    
    if ligand is None or receptor is None:
        raise ValueError("Both 'ligand' and 'receptor' must be specified.")
    
    pair = f"{ligand}{lr_sep}{receptor}"
    mask = adata.var_names.str.endswith(pair)
    if not mask.any():
        raise ValueError(f"❌ No LR features found for: {pair} (e.g., 'Sender{lr_sep}{pair}'). Check ligand/receptor names or 'lr_sep'.")

    # Use a view or copy for safety and performance
    sub = adata[:, mask].copy()

    # 2.Parsing of Variable Names
    # Expects format: 'sender^ligand^receptor'
    try:
        parts = sub.var_names.to_series().str.rsplit(lr_sep, n=2, expand=True)
        if parts.shape[1] != 3:
             # Fallback/Error if the split didn't yield 3 parts
            raise ValueError(f"Could not parse var_names using separator '{lr_sep}'. Expected format: 'sender{lr_sep}ligand{lr_sep}receptor'.") from None

        parts.columns = ["sender", "ligand", "receptor"]
        # Assign the parsed sender labels to `sub.var`
        sub.var["sender"] = parts["sender"].values
    except ValueError as e:
        raise ValueError(f"Error during var_names parsing: {e}") from e

    # 3. Prepare Data for Heatmap
    receiver = sub.obs[groupby].astype(str) # Convert to string for consistent grouping

    # Extract expression data
    X = sub.X.toarray() if hasattr(sub.X, "toarray") else sub.X
    df = pd.DataFrame(X, index=sub.obs_names, columns=sub.var["sender"])
    df["receiver"] = receiver.values

    # Calculate mean inflow (mean signal for each sender→receiver interaction)
    # The result is (sender rows) x (receiver columns)
    heatmap_df = df.groupby("receiver").mean().T

    if filter_fun is not None:
        heatmap_df = filter_fun(heatmap_df)
        if heatmap_df.empty:
            raise ValueError("⚠️ After filter_fun, the heatmap is empty.")

    # 5. Ordering for Visual Clarity
    mat = heatmap_df.copy()

    # Order rows (senders) and columns (receivers) by their total signal sum
    row_order = mat.sum(axis=1).sort_values(ascending=False).index
    col_order = mat.sum(axis=0).sort_values(ascending=False).index
    mat = mat.loc[row_order, col_order]

    # 6. Clipping (Color Scale Stabilization)
    vmin = vmax = None
    if clip:
        # Calculate 5th and 95th percentiles from non-NaN values
        valid_data = mat.values[~np.isnan(mat.values)]
        if valid_data.size > 0:
            vmin = np.nanquantile(valid_data, 0.05)
            vmax = np.nanquantile(valid_data, 0.95)
        else:
            print("Warning: Data is empty or all NaN, skipping clipping.")
            clip = False 

    # 7. Plotting
    plt.figure(figsize=figure_size)

    # Set default seaborn heatmap aesthetics for a clean look
    default_kwargs = {
        'linewidths': 0.4,
        'linecolor': "white",
        'square': True,
        'cbar_kws': {'shrink': 0.8} # Shrink the colorbar slightly
    }

    # Merge default kwargs with user-provided kwargs, prioritizing user input
    plot_kwargs = {**default_kwargs, **kwargs}

    sns.heatmap(
        mat,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        **plot_kwargs
    )

    # 8. Final Touches
    plt.xlabel(xlabel, fontsize=12, fontweight='bold')
    plt.ylabel(ylabel, fontsize=12, fontweight='bold')
    plt.title(title or f"Sender→Receiver Average Inflow Score for {ligand}-{receptor}", fontsize=14, fontweight='bold')
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.show()

    return mat
