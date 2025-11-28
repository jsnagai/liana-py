from __future__ import annotations
import anndata
import pandas
from plotnine import (
    aes,
    element_rect,
    element_text,
    facet_grid,
    geom_point,
    ggplot,
    labs,
    scale_color_cmap,
    scale_size_continuous,
    theme,
    theme_bw,
)

import seaborn as sns
import matplotlib.pyplot as plt
from liana._constants import DefaultValues as V
from liana._constants import Keys as K
from liana._docs import d
from liana.plotting._common import _check_var, _filter_by, _get_top_n, _invert_scores, _prep_liana_res

@d.dedent
def heatmap(
    adata: anndata.AnnData = None,
    uns_key: str = K.uns_key,
    liana_res: pandas.DataFrame = None,
    ligand: str = None,
    receptor: str = None,
    groupby:  str = None,
    pval_annotation: str = None,
    cmap: str = V.cmap,
    figure_size: tuple = (10, 10),
    filter_fun: callable = None,
):
    
    """
    Heatmap of ligand-receptor (LR) interactions by source and target cells

    Parameters
    ----------
    %(adata)s
    %(uns_key)s
    %(liana_res)s
    ligand : str
        Name of the ligand gene.
    receptor : str
        Name of the receptor gene.
    %(groupby)s
    pval_annotation : str | None, optional
        Type of p-value annotation to display on the heatmap. Options are 'star' for
        asterisks, 'number' for numerical values, and 'none' for no annotation. Defaults to None.
    %(cmap)s
    %(filter_fun)s
    %(figure_size)s

    """
    
    if not isinstance(adata, anndata.AnnData):
        raise TypeError("adata must be an AnnData object.")
    if groupby not in adata.obs.columns:
        raise ValueError(f"'{groupby}' not found in adata.obs. Available columns are: {list(adata.obs.columns)}")
    if ligand is None or receptor is None:
        raise ValueError("Both 'ligand' and 'receptor' must be specified.")
    if uns_key not in adata.uns:
        raise ValueError(f"'{uns_key}' not found in adata.uns. Available keys are: {list(adata.uns.keys())}")
    

    #Consider all unique labels in the groupby column as source and target labels
    source_labels = adata.obs[groupby].unique().tolist()
    target_labels = adata.obs[groupby].unique().tolist()

    liana_res = _prep_liana_res(adata=adata,
                            liana_res=liana_res,
                            source_labels=source_labels,
                            target_labels=target_labels,
                            ligand_complex=ligand,
                            receptor_complex=receptor,
                            uns_key=uns_key)
    
    
    liana_res = _filter_by(liana_res, filter_fun)
    
    df = liana_res[(liana_res["ligand"] == ligand) & (liana_res["receptor"] == receptor)].copy()
    if df.empty:
         raise ValueError(f"No rows found for ligand={ligand}, receptor={receptor} in liana_res")


    if filter_fun is not None:
        # 1. Identify which rows satisfy the filter function
        mask_satisfies_filter = df.apply(filter_fun, axis=1).astype(bool)
        df_satisfy = df[mask_satisfies_filter]
        if df_satisfy.empty:
            raise ValueError(f"No interactions found for ligand={ligand}, receptor={receptor} after filter_fun is applied.")

        # 2. Collect the 'source' and 'target' labels that satisfied the filter
        # Keep all rows where the source AND the target has at least one filtered interaction
        sources_to_keep = df_satisfy["source"].unique()
        targets_to_keep = df_satisfy["target"].unique()
        relevant_cell_types = set(sources_to_keep) | set(targets_to_keep)
        # 3. Apply the new filtering rule to the original dataframe
        rows_to_keep_mask = (df["source"].isin(relevant_cell_types)) & (df["target"].isin(relevant_cell_types))
        df = df[rows_to_keep_mask]
        if df.empty:
            # This should technically not be hit if df_satisfy wasn't empty, but kept for robustness
            raise ValueError(f"No interactions found for ligand={ligand}, receptor={receptor} after groupby-level filtering.")

    heatmap_df = df.pivot(index="source", columns="target", values="lr_mean")

    # convert p-values to stars
    def p_to_star(p):
        if p < 0.001: return "***"
        if p < 0.01:  return "**"
        if p < 0.05:  return "*"
        return ""

    # pval annotation 
    if pval_annotation is None or pval_annotation == "none":
        df["annot"] = ""
    elif pval_annotation == "star":
        df["annot"] = df["pval"].apply(p_to_star)
    elif pval_annotation == "number":
        df["annot"] = df["pval"].apply(lambda x: f"{x:.3g}")

    pval_m = df.pivot(index="source", columns="target", values="annot")


    plt.figure(figsize=figure_size)
    sns.heatmap(
        heatmap_df,
        annot=pval_m,
        fmt="s",
        cmap=cmap,
        cbar_kws={"label": "lr_mean"},
        linewidths=0.5,
        linecolor="gray"
    )

    plt.title(df["interaction"].iloc[0])
    plt.xlabel("Target cell type")
    plt.ylabel("Source cell type")
    plt.show()