import numpy as np
import pytest
from scipy.sparse import csr_matrix

from liana.method import inflow
from liana.testing._sample_anndata import generate_toy_spatial
from liana.testing._sample_resource import sample_resource
from liana.utils.transform import zi_minmax


@pytest.fixture
def spatial_adata():
    """Generate spatial AnnData with connectivity matrix."""
    adata = generate_toy_spatial()
    # Ensure spatial_connectivities exists
    if 'spatial_connectivities' not in adata.obsp:
        from liana.utils import spatial_neighbors
        spatial_neighbors(adata, bandwidth=100, spatial_key='spatial')
    return adata


def test_inflow_basic_structure(spatial_adata):
    """Test basic execution and output structure preservation."""
    lrdata = inflow(
        spatial_adata,
        groupby='bulk_labels',
        resource_name='consensus',
        use_raw=True
    )

    # Check output structure
    assert isinstance(lrdata, type(spatial_adata))
    assert lrdata.shape == (spatial_adata.shape[0], 323)  # Fixed expected shape
    assert lrdata.obs.shape == spatial_adata.obs.shape
    assert lrdata.obsm.keys() == spatial_adata.obsm.keys()

    # Check var index format: "celltype^ligand^receptor"
    assert all('^' in idx for idx in lrdata.var_names)
    for idx in lrdata.var_names:
        parts = idx.split('^')
        assert len(parts) == 3
        celltype, ligand, receptor = parts
        assert celltype in spatial_adata.obs['bulk_labels'].unique()

    # Check sparse matrix
    assert isinstance(lrdata.X, csr_matrix)

    # Check obs and obsm preserved
    assert lrdata.obs.equals(spatial_adata.obs)
    assert 'spatial' in lrdata.obsm
    np.testing.assert_array_equal(
        lrdata.obsm['spatial'],
        spatial_adata.obsm['spatial']
    )

    # Check obsp preserved
    assert 'spatial_connectivities' in lrdata.obsp


def test_inflow_with_transform(spatial_adata):
    """Test inflow with zi_minmax transformation."""
    lrdata = inflow(
        spatial_adata,
        groupby='bulk_labels',
        resource_name='consensus',
        x_transform=zi_minmax,
        y_transform=zi_minmax,
        use_raw=False
    )

    assert lrdata.shape[0] == spatial_adata.shape[0]
    # Transformed values should be in [0, 1]
    assert lrdata.X.min() >= 0
    assert lrdata.X.max() <= 1


def test_inflow_nz_prop_filter(spatial_adata):
    """Test filtering by non-zero proportion."""
    # Strict filter
    lrdata_strict = inflow(
        spatial_adata,
        groupby='bulk_labels',
        resource_name='consensus',
        nz_prop=0.2,
        use_raw=True
    )

    # Lenient filter
    lrdata_lenient = inflow(
        spatial_adata,
        groupby='bulk_labels',
        resource_name='consensus',
        nz_prop=0.001,
        use_raw=True
    )

    # Strict filter should have fewer or equal interactions
    assert lrdata_strict.shape[1] <= lrdata_lenient.shape[1]


def test_inflow_custom_resource(spatial_adata):
    """Test with custom resource DataFrame."""
    resource = sample_resource(spatial_adata, n_lrs=10)

    lrdata = inflow(
        spatial_adata,
        groupby='bulk_labels',
        resource=resource,
        use_raw=True
    )

    assert lrdata.shape[1] > 0
    assert lrdata.shape[0] == spatial_adata.shape[0]


def test_inflow_numerical_correctness(spatial_adata):
    """Test numerical correctness of inflow scores."""
    lrdata = inflow(
        spatial_adata,
        groupby='bulk_labels',
        resource_name='consensus',
        use_raw=True
    )

    # Check specific numerical values (regression test)
    np.testing.assert_almost_equal(lrdata.X.mean(), 0.041507, decimal=3)  # Replace with actual
    np.testing.assert_almost_equal(lrdata.X.sum(), 9384.73809, decimal=3)    # Replace with actual


def test_inflow_missing_connectivity(spatial_adata):
    """Test error when spatial_connectivities is missing."""
    # Remove spatial connectivity
    del spatial_adata.obsp['spatial_connectivities']

    with pytest.raises(ValueError, match="spatial_connectivities"):
        inflow(
            spatial_adata,
            groupby='bulk_labels',
            resource_name='consensus',
            use_raw=True
        )


def test_inflow_no_features_pass_filter(spatial_adata):
    """Test error when no features pass nz_prop filter."""
    with pytest.raises(ValueError, match="No features passed"):
        inflow(
            spatial_adata,
            groupby='bulk_labels',
            resource_name='consensus',
            nz_prop=0.99,  # Very strict filter
            use_raw=True
        )


def test_inflow_invalid_groupby(spatial_adata):
    """Test error with invalid groupby column."""
    with pytest.raises(KeyError):
        inflow(
            spatial_adata,
            groupby='nonexistent_column',
            resource_name='consensus',
            use_raw=True
        )
