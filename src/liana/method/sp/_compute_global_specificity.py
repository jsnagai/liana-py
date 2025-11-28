import numpy as np
import pandas as pd
from anndata import AnnData
from joblib import Parallel, delayed
from scipy.sparse import csr_matrix, diags, issparse
from tqdm import trange, tqdm

from liana._constants import DefaultValues as V
from liana._logging import _logg
from liana._docs import d


def _split_complex(name: str, complex_sep: str = V.complex_sep):
    """Helper for splitting complex names."""
    toks = name.split(complex_sep)
    if len(toks) > 1:
        return toks[0], complex_sep.join(toks[1:])
    else:
        return name, name

def _run_single_permutation(
    shuffled_labels: np.ndarray,
    X,
    celltype_names: list,
    interaction_names: np.ndarray,
    lr_sep: str = V.lr_sep,
) -> dict:
    """Run a single permutation using precomputed shuffled labels."""
    cats = pd.Categorical(shuffled_labels, categories=celltype_names, ordered=False)
    col_idx = cats.codes

    n_cells = col_idx.size
    n_types = len(celltype_names)
    row_idx = np.arange(n_cells)

    # Build sparse one-hot and normalize columns
    t_sparse = csr_matrix(
        (np.ones(n_cells, dtype=np.float64), (row_idx, col_idx)),
        shape=(n_cells, n_types)
    )
    col_sums = np.asarray(t_sparse.sum(axis=0)).ravel()
    col_sums[col_sums == 0] = 1.0
    t_sparse = t_sparse @ diags(1.0 / col_sums)

    # Ensure X is sparse CSR
    X = X if issparse(X) else csr_matrix(X)

    # Aggregation
    result = t_sparse.T @ X
    values = (result.toarray() if issparse(result) else np.asarray(result)).ravel()

    # Build keys with vectorized string ops
    interaction_names_tiled = np.tile(interaction_names.astype(str), n_types)
    celltype_names_repeated = np.repeat(
        np.asarray(celltype_names, dtype=str),
        interaction_names.shape[0]
    )
    keys = np.char.add(np.char.add(interaction_names_tiled, lr_sep), celltype_names_repeated)

    return dict(zip(keys.tolist(), values.tolist(), strict=True))

@d.dedent
def compute_global_specificity(
    adata: AnnData,
    groupby: str,
    lr_sep: str = V.lr_sep,
    complex_sep: str = V.complex_sep,
    n_perms: int = V.n_perms,
    seed: int = V.seed,
    n_jobs: int = -1,
    verbose: bool = V.verbose,
) -> None:
    """
    Computes global specificity and calculates permutation test p-values.

    Args:
        %(adata)s
        %(groupby)s
        lr_sep (str, optional): Separator for names. Defaults to `V.lr_sep` ('^').
        complex_sep (str, optional): Separator for splitting complex names. Defaults to "_".
        %(n_perms)s
        %(seed)s
        %(verbose)s

    Returns
    -------
        None: The result with 'lr_mean' and 'pval' is stored in `adata.uns["global_interactions"]`.
    """
    if groupby not in adata.obs.columns:
        raise KeyError(
            f"`groupby`='{groupby}' not found in adata.obs. "
            "Use the same grouping column used to build adata."
        )

    rng_main = np.random.default_rng(seed)
    original_groupby_labels = adata.obs[groupby].copy()
    # --- Part A: Observed score (sparse) ---

    # Fixed column order for cell types
    celltypes = pd.get_dummies(adata.obs[groupby])
    celltype_names = list(celltypes.columns)

    # Map original labels to fixed indices
    cats_obs = pd.Categorical(adata.obs[groupby].values, categories=celltype_names, ordered=False)
    col_idx_obs = cats_obs.codes

    n_cells = col_idx_obs.size
    n_types = len(celltype_names)
    row_idx = np.arange(n_cells)

    # Sparse one-hot and normalize columns
    t_sparse = csr_matrix(
        (np.ones(n_cells, dtype=np.float64), (row_idx, col_idx_obs)),
        shape=(n_cells, n_types)
    )
    col_sums = np.asarray(t_sparse.sum(axis=0)).ravel()
    col_sums[col_sums == 0] = 1.0
    t_sparse = t_sparse @ diags(1.0 / col_sums)

    # Ensure X is sparse CSR
    X_raw = adata.X
    X = X_raw if issparse(X_raw) else csr_matrix(X_raw)

    # Aggregation
    result = t_sparse.T @ X
    values = (result.toarray() if issparse(result) else np.asarray(result)).ravel()

    # Names
    interaction_names = adata.var.index.astype(str).values

    # Build full names "L^R^Source^Target"
    full_names = np.char.add(
        np.char.add(
            np.tile(interaction_names.astype(str), n_types),
            lr_sep
        ),
        np.repeat(np.asarray(celltype_names, dtype=str), interaction_names.shape[0])
    )

    # Build observed df
    parts = [n.split(lr_sep) for n in full_names.tolist()]
    df = pd.DataFrame(parts, columns=["source", "ligand", "receptor", "target"])
    df["lr_mean"] = values
    df["pval"] = np.nan

    # Complex parsing
    lig_primary, lig_complex = zip(*df["ligand"].map(lambda x: _split_complex(x, complex_sep)), strict=True)
    rec_primary, rec_complex = zip(*df["receptor"].map(lambda x: _split_complex(x, complex_sep)), strict=True)
    df["ligand"] = lig_primary
    df["ligand_complex"] = lig_complex
    df["receptor"] = rec_primary
    df["receptor_complex"] = rec_complex

    observed_df = df.copy()

    # Prepare for permutation test
    interaction_keys_for_init = full_names.tolist()
    observed_score_map = dict(zip(interaction_keys_for_init, observed_df["lr_mean"].values, strict=True))
    perm_matrix = {key: [] for key in interaction_keys_for_init}

    # --- Part B: Permutation test ---

    # Precompute all permutations
    permuted_labels_list = [
        rng_main.permutation(original_groupby_labels.values)
        for _ in range(n_perms)
    ]

    if verbose:
        _logg(f"Running {n_perms} permutations in parallel...")

    joblib_verbose = 5 if verbose else 0

    # Run in parallel
    permuted_scores_list = Parallel(n_jobs=n_jobs, verbose=joblib_verbose)(
        delayed(_run_single_permutation)(
            shuffled_labels=shuffled_labels,
            X=X,
            interaction_names=interaction_names,
            celltype_names=celltype_names,
            lr_sep=lr_sep,
        ) for shuffled_labels in tqdm(permuted_labels_list, desc="Running permutations")
    )

    # Aggregate
    for scores_dict in permuted_scores_list:
        for key, score in scores_dict.items():
            perm_matrix[key].append(score)

    # Compute p-values
    n = float(n_perms + 1)
    pvals = []
    for key in interaction_keys_for_init:
        real_score = observed_score_map.get(key)
        perm_scores = np.asarray(perm_matrix[key], dtype=np.float64)
        pval = (np.sum(perm_scores >= real_score) + 1.0) / n
        pvals.append(pval)

    observed_df["pval"] = pvals
    observed_df = observed_df[[
        "ligand", "ligand_complex", "receptor", "receptor_complex",
        "source", "target", "lr_mean", "pval"
    ]]

    # Save result
    adata.uns["global_interactions"] = observed_df
