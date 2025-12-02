import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from numpy.testing import assert_allclose
from liana.method.sp._compute_global_specificity import compute_global_specificity

def _make_example_adata():
    # 4 cells, 2 interactions
    # Var names are structured so that splitting by lr_sep="^" yields 3 parts in the interaction name,
    # and adding the celltype produces the 4th part expected by the function.
    var_names = ["SRC1^L1^R1", "SRC2^L2^R2"]
    X = np.array([
        [1.0, 2.0],  # cell 0, group A
        [3.0, 4.0],  # cell 1, group A
        [5.0, 6.0],  # cell 2, group B
        [7.0, 8.0],  # cell 3, group B
    ], dtype=float)
    obs = pd.DataFrame({"ct": ["A", "A", "B", "B"]})
    var = pd.DataFrame(index=var_names)
    adata = AnnData(X=X, obs=obs, var=var)
    return adata


def test_compute_global_specificity_basic():
    adata = _make_example_adata()
    # Run with small number of permutations and single job for determinism in tests
    compute_global_specificity(
        adata,
        groupby="ct",
        lr_sep="^",
        complex_sep="_",
        n_perms=3,
        seed=0,
        n_jobs=1,
        verbose=False,
    )

    assert "global_interactions" in adata.uns
    df = adata.uns["global_interactions"]
    # Expected columns order
    expected_cols = [
        "ligand", "ligand_complex", "receptor", "receptor_complex",
        "source", "target", "specificity_rank", "pval"
    ]
    assert list(df.columns) == expected_cols

    # There are 2 interactions * 2 cell types = 4 rows
    assert df.shape[0] == 4

    # Specificity ranks should equal the group-wise means:
    # For group A (cells 0,1): interaction1 -> (1+3)/2=2, interaction2 -> (2+4)/2=3
    # For group B (cells 2,3): interaction1 -> (5+7)/2=6, interaction2 -> (6+8)/2=7
    expected_specificity = [2.0, 3.0, 6.0, 7.0]
    assert_allclose(df["specificity_rank"].values.astype(float), expected_specificity)

    # p-values must be between 0 and 1
    pvals = df["pval"].values.astype(float)
    assert np.all(pvals >= 0.0) and np.all(pvals <= 1.0)

    # Check that ligand/receptor parsing worked (from "SRC^L^R^Target" -> ligand == L, receptor == R)
    # Because var names were "SRC1^L1^R1" and celltypes "A"/"B", ligand should be "L1"/"L2" etc.
    assert set(df["ligand"].unique()) == {"L1", "L2"}
    assert set(df["receptor"].unique()) == {"R1", "R2"}


def test_compute_global_specificity_invalid_groupby_raises():
    adata = _make_example_adata()
    with pytest.raises(KeyError):
        compute_global_specificity(adata, groupby="nonexistent", n_perms=1, n_jobs=1, verbose=False)