import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import trim_mean
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize

from liana._constants import DefaultValues as V
from liana._constants import Keys as K
from liana._docs import d


def _gaussian(distance_mtx, bandwidth):
    return np.exp(-(distance_mtx ** 2.0) / (2.0 * bandwidth ** 2.0))

def _misty_rbf(distance_mtx, bandwidth):
    return np.exp(-(distance_mtx ** 2.0) / (bandwidth ** 2.0))

def _exponential(distance_mtx, bandwidth):
    return np.exp(-distance_mtx / bandwidth)

def _linear(distance_mtx, bandwidth):
    connectivity = 1 - distance_mtx / bandwidth
    return np.clip(connectivity, a_min=0, a_max=np.inf)


def _kernel_function(distance_mtx, bandwidth, kernel):
    families = ['gaussian', 'exponential', 'linear', 'misty_rbf']
    if kernel not in families:
        raise AssertionError(f"{kernel} must be a member of {families}")

    if kernel == 'gaussian':
        return _gaussian(distance_mtx, bandwidth)
    elif kernel == 'misty_rbf':
        return _misty_rbf(distance_mtx, bandwidth)
    elif kernel == 'exponential':
        return _exponential(distance_mtx, bandwidth)
    elif kernel == 'linear':
        return _linear(distance_mtx, bandwidth)
    else:
        raise ValueError("Please specify a valid family to generate connectivity weights")


@d.dedent
def spatial_neighbors(adata: AnnData,
                      bandwidth=None,
                      cutoff=0.1,
                      max_neighbours=100,
                      kernel='gaussian',
                      set_diag=False,
                      zoi=0,
                      standardize=False,
                      reference=None,
                      spatial_key=K.spatial_key,
                      key_added=K.spatial_key,
                      inplace=V.inplace
                      ):
    """
    Generate spatial connectivity weights using Euclidean distance.

    Parameters
    ----------
    %(adata)s
    bandwidth
         Denotes signaling length (`l`) and controls the maximum distance at which two spots are considered.
         Corresponds to the units in which spatial coordinates are expressed.
    cutoff
        Values below this cutoff will be set to 0.
    max_neighbours
        Maximum nearest neighbours to be considered when generating spatial connectivity weights.
        Essentially, the maximum number of edges in the spatial connectivity graph.
    %(kernel)s
    set_diag
        Logical, sets connectivity diagonal to 0 if `False`. Default is `True`.
    zoi
        Zone of indifference. Values below this cutoff will be set to `np.inf`.
    standardize
        Whether to (l1) standardize spatial proximities (connectivities) so that they sum to 1.
        This plays a role when weighing border regions prior to downstream methods, as the number of spots
        in the border region (and hence the sum of proximities) is smaller than the number of spots in the center.
        Relevant for methods with unstandardized scores (e.g. product). Default is `False`.
    reference
        Reference coordinates to use when generating spatial connectivity weights.
        If `None`, uses the spatial coordinates in `adata.obsm[spatial_key]`.
        This is only relevant if you want to use a different set of coordinates to generate spatial connectivity weights.
    %(spatial_key)s
    key_added
        Key to add to `adata.obsp` if `inplace = True`. If reference is not `None`, key will be added to `adata.obsm`.
    %(inplace)s

    Notes
    -----
    This function is adapted from mistyR, and is set to be consistent with
    the `squidpy.gr.spatial_neighbors` function in the `squidpy` package.

    Returns
    -------
    If ``inplace = False``, returns an `np.array` with spatial connectivity weights.
    Otherwise, modifies the ``adata`` object with the following key:
        - :attr:`anndata.AnnData.obsp` ``['{key_added}_connectivities']`` with the aforementioned array

    """
    if cutoff is None:
        raise ValueError("`cutoff` must be provided!")
    assert spatial_key in adata.obsm
    families = ['gaussian', 'exponential', 'linear', 'misty_rbf']
    if kernel not in families:
        raise AssertionError(f"{kernel} must be a member of {families}")
    if bandwidth is None:
        raise ValueError("Please specify a bandwidth")

    coordinates = adata.obsm[spatial_key]

    if reference is None:
        _reference = coordinates
    else:
        _reference = reference

    tree = NearestNeighbors(n_neighbors=max_neighbours + 1, # +1 to exclude self
                            algorithm='ball_tree',
                            metric='euclidean').fit(_reference)
    dist = tree.kneighbors_graph(coordinates, mode='distance')

    # prevent float overflow
    bandwidth = np.array(bandwidth, dtype=np.float64)

    # define zone of indifference
    dist.data[dist.data < zoi] = np.inf

    # NOTE: dist gets converted to a connectivity (proximity) matrix
    dist.data = _kernel_function(dist.data, bandwidth, kernel)

    if not set_diag:
        dist.setdiag(0)
    if cutoff is not None:
        dist.data = dist.data * (dist.data > cutoff)
    if standardize:
        dist = normalize(dist, axis=1, norm='l1')

    spot_n = dist.shape[0]
    if reference is None:
        assert spot_n == adata.shape[0]
    if spot_n > 1000:
        dist = dist.astype(np.float32)

    if inplace:
        if reference is not None:
            adata.obsm[f'{key_added}_connectivities'] = dist
        else:
            adata.obsp[f'{key_added}_connectivities'] = dist

    return None if inplace else dist


@d.dedent
def cell_type_spatial_proximity(
    groupby,
    coordinates,
    interaction_range=250,
    contact_range=None,
    min_cells_in_proximity=10,
    trim_fraction=0.1,
    kernel='gaussian'
):
    """
    Computes aggregated spatial statistics and proximity scores between cell types.

    This function calculates pairwise proximity between cell types based on nearest neighbor
    distances in spatial coordinates. It returns a DataFrame with proximity scores that can
    be used to weight ligand-receptor interactions by spatial co-localization.

    Parameters
    ----------
    %(groupby)s
    %(coordinates)s
    interaction_range : float, optional
        Bandwidth for 'proximity' score calculation and distance threshold for 'interacting' flag.
        Default is 250.
    contact_range : float, optional
        Bandwidth for 'contact_proximity' calculation and distance threshold for 'contact_interacting' flag.
        If None, contact proximity is not calculated. Default is None.
    min_cells_in_proximity : int, optional
        Minimum number of cell pairs within range required to flag an interaction as significant.
        Default is 10.
    trim_fraction : float, optional
        Fraction of outliers to trim from each tail when calculating mean distance (0-0.5).
        Default is 0.1 (trim 10% from each tail).
    %(kernel)s

    Returns
    -------
    pd.DataFrame
        DataFrame with columns:
        - source: source cell type
        - target: target cell type
        - mean_distance: trimmed mean distance between cell types
        - interacting: binary flag (1 if >= min_cells_in_proximity pairs within interaction_range, else 0)
        - proximity: proximity score calculated by applying kernel to mean_distance with interaction_range bandwidth
        - contact_interacting: (optional, if contact_range is not None) binary flag for contact-range interactions
        - contact_proximity: (optional, if contact_range is not None) proximity score using contact_range bandwidth

    Notes
    -----
    - Performance scales as O(n_cell_types² × n_cells), which is acceptable for typical datasets
      (5-30 cell types) but may be slower with 100+ cell types.
    - Self-interactions exclude the cell itself as its own neighbor to avoid zero distances.
    - Missing proximity values (e.g., cell types that never co-localize) will result in NaN,
      which should be filled with 0.0 when merging with interaction results.

    Examples
    --------
    >>> import numpy as np
    >>> import pandas as pd
    >>> # Create example data
    >>> groupby = np.array(['TypeA', 'TypeA', 'TypeB', 'TypeB'])
    >>> coordinates = np.array([[0, 0], [1, 1], [10, 10], [11, 11]])
    >>> proximity_df = cell_type_spatial_proximity(groupby, coordinates)
    >>> proximity_df
    """
    groupby = np.asarray(groupby)
    coordinates = np.asarray(coordinates, dtype=float)

    unique_types = np.unique(groupby)
    stats_list = []

    # Iterate through all cell type pairs
    for type_a in unique_types:
        idx_a = np.where(groupby == type_a)[0]
        coords_a = coordinates[idx_a]

        for type_b in unique_types:
            idx_b = np.where(groupby == type_b)[0]
            coords_b = coordinates[idx_b]

            if len(idx_a) == 0 or len(idx_b) == 0:
                continue

            # Handle self-interaction (exclude cell itself as neighbor)
            is_self = (type_a == type_b)
            k_neighbors = 2 if is_self else 1

            if is_self and len(idx_b) < 2:
                continue

            # Nearest neighbor search (1-NN)
            nn = NearestNeighbors(n_neighbors=k_neighbors, metric="euclidean", n_jobs=-1)
            nn.fit(coords_b)
            distances, _ = nn.kneighbors(coords_a)

            # If self, take 2nd column; if different, take 1st column
            raw_dists = distances[:, 1] if is_self else distances[:, 0]

            # --- Aggregation ---

            # 1. Trimmed mean distance (core metric)
            avg_dist = trim_mean(raw_dists, proportiontocut=trim_fraction)

            # 2. Binary flags (significance based on counts)
            count_long = np.sum(raw_dists <= interaction_range)
            is_interacting = count_long >= min_cells_in_proximity

            # 3. Proximity score (kernel applied to mean_distance)
            prox_score = _kernel_function(avg_dist, bandwidth=interaction_range, kernel=kernel)

            # Build result dict
            result_dict = {
                "source": type_a,
                "target": type_b,
                "mean_distance": avg_dist,
                "interacting": int(is_interacting),
                "proximity": prox_score
            }

            # 4. Optional contact proximity
            if contact_range is not None:
                count_short = np.sum(raw_dists <= contact_range)
                is_physically_interacting = count_short >= min_cells_in_proximity
                contact_prox_score = _kernel_function(avg_dist, bandwidth=contact_range, kernel=kernel)

                result_dict["contact_interacting"] = int(is_physically_interacting)
                result_dict["contact_proximity"] = contact_prox_score

            stats_list.append(result_dict)

    return pd.DataFrame(stats_list)
