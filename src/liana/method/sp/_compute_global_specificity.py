import numpy as np
import pandas as pd
from anndata import AnnData
from joblib import Parallel, delayed
from scipy.sparse import csr_matrix, diags, issparse
from tqdm import trange, tqdm

from liana._constants import DefaultValues as V
from liana._logging import _logg
from liana._docs import d
from liana.method._pipe_utils._pre import _choose_mtx_rep

#Function to compute group means
def get_group_mean(X, groupby_labels, var_names, groups_order=None):
    s = pd.Series(groupby_labels)
    groups_dum = pd.get_dummies(s, dummy_na=False)
    if groups_order is not None:
        groups_dum = groups_dum.reindex(columns=groups_order, fill_value=0)
    t_sparse = csr_matrix(groups_dum)
    col_sums = np.asarray(t_sparse.sum(axis=0)).ravel()
    col_sums[col_sums == 0] = 1.0
    t_sparse = t_sparse @ diags(1.0 / col_sums)
    result = t_sparse.T @ X
    arr = result.toarray()
    df_raw = pd.DataFrame(arr, index=groups_dum.columns, columns=var_names)
    return df_raw


@d.dedent
def compute_global_specificity(
    adata: AnnData,
    groupby: str,
    lr_sep: str | None = V.lr_sep,
    n_perms: int = V.n_perms,
    seed: int = V.seed,
    n_jobs: int = -1,
    verbose: bool = V.verbose,
    use_raw: bool = V.use_raw,
    layer: str = V.layer,
    uns_key: str = "global_interactions",
) -> None:
    """
    Computes global specificity and calculates permutation test p-values.

    Args:
        %(adata)s
        %(groupby)s
        %(lr_sep)s
        %(n_perms)s
        %(seed)s
        n_jobs (int, optional): Number of parallel jobs. Defaults to -1 (all available cores).
        %(verbose)s
        %(use_raw)s
        %(layer)s
        %(uns_key)s

    Returns
    -------
        None: The result with 'lr_mean' and 'pval' is stored in `adata.uns["global_interactions"]`.
    """
    if groupby not in adata.obs.columns:
        raise KeyError(f"`groupby`='{groupby}' not found in adata.obs.")
    
    X = _choose_mtx_rep(adata, layer=layer, use_raw=use_raw)
    var_names = adata.var_names
    original_groupby_labels = adata.obs[groupby].astype('category')
    groups_order = list(original_groupby_labels.cat.categories)

    #Compute observed statistic
    df_obs = get_group_mean(X, original_groupby_labels, var_names, groups_order=groups_order)
    obs_values = df_obs.values

    if verbose:
        _logg(f"Running {n_perms} permutations in parallel...")

    joblib_verbose = 5 if verbose else 0

    rng_main = np.random.default_rng(seed)
    perm_stats_list = Parallel(n_jobs=n_jobs, verbose=joblib_verbose)(
        delayed(get_group_mean)(X, rng_main.permutation(original_groupby_labels.values), var_names, groups_order=groups_order)
        for _ in tqdm(range(n_perms), desc="Running permutations"))

    #Convert results to array
    perm_stats = np.array([df_perm.values for df_perm in perm_stats_list])

    # Compute empirical p-values
    k = np.sum(np.abs(perm_stats) >= np.abs(obs_values), axis=0)
    p_values = (k + 1) / (n_perms + 1)
    pval_df = pd.DataFrame(p_values, index=df_obs.index, columns=df_obs.columns)
    pval_melted = pval_df.reset_index().melt(id_vars='index', var_name='feature', value_name='pval')
    
    #Final df
    df_melted = df_obs.reset_index().melt(id_vars='index', var_name='feature', value_name='lr_mean')
    merged_df = df_melted.merge(pval_melted[['index', 'feature', 'pval']], on=['index', 'feature'], how='left')
    
    if lr_sep is not None:
        merged_df[['source', 'ligand_complex', 'receptor_complex']] = merged_df['feature'].str.split(lr_sep, expand=True)
        merged_df = merged_df.rename(columns={'index': 'target'})
        final_df = merged_df[['source', 'ligand_complex', 'receptor_complex', 'target', 'lr_mean', 'pval']]
    else:
        final_df = merged_df

    final_df = final_df.sort_values('pval', ascending=True)  

    #Save result
    adata.uns[uns_key] = final_df