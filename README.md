# LIANA+: an all-in-one cell-cell communication framework <img src="https://raw.githubusercontent.com/saezlab/liana-py/main/docs/_static/logo.png" align="right" height="125">

<!-- badges: start -->
[![main](https://github.com/saezlab/liana-py/actions/workflows/main.yml/badge.svg)](https://github.com/saezlab/liana-py/actions)
[![GitHub issues](https://img.shields.io/github/issues/saezlab/liana-py.svg)](https://github.com/saezlab/liana-py/issues/)
[![Documentation Status](https://readthedocs.org/projects/liana-py/badge/?version=latest)](https://liana-py.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/saezlab/liana-py/branch/main/graph/badge.svg?token=TM0P29KKN5)](https://codecov.io/gh/saezlab/liana-py)
[![Downloads](https://static.pepy.tech/badge/liana)](https://pepy.tech/project/liana)
<!-- badges: end -->

LIANA+ is a scalable framework that adapts and extends existing methods and knowledge to study cell-cell communication in single-cell, spatially-resolved, and multi-modal omics data. It is part of the [scverse ecosystem](https://github.com/scverse), and relies on [AnnData](https://github.com/scverse/anndata) & [MuData](https://github.com/scverse/mudata) objects as input.

<img src="https://raw.githubusercontent.com/saezlab/liana-py/main/docs/_static/abstract.png" width="700" align="center">

## Contributions

We welcome suggestions, ideas, and contributions! Please do not hesitate to contact us, open issues, and check the [contributions guide](https://liana-py.readthedocs.io/en/latest/contributing.html).

## Vignettes
A set of extensive vignettes can be found in the [LIANA+ documentation](https://liana-py.readthedocs.io/en/latest/).

## Decision Tree

```mermaid
flowchart TD
    Start[What type of data?] --> Spatial{Spatial<br/>coordinates?}
    Start --> Modal{Multi-modal?}

    %% Spatial branch
    Spatial -->|Yes| SpatialRes{Resolution?}
    SpatialRes -->|Single-cell| Inflow[Inflow Score]
    SpatialRes -->|Spot-based| SpatialType{Analysis type?}
    SpatialType -->|Bivariate| LocalQ{Local<br/>interactions?}
    LocalQ -->|Yes| Local[Local Bivariate Metrics]
    LocalQ -->|No| Global[Global Bivariate Metrics]
    SpatialType -->|Unsupervised| MISTy[Multi-view Learning]

    %% Non-spatial branch
    Spatial -->|No| Compare{Compare across<br/>samples?}
    Compare -->|Yes| Contrast{Specific<br/>contrast?}
    Contrast -->|Yes| Targeted[Differential Contrasts]
    Contrast -->|No| MOFA[MOFA+]
    Contrast -->|No| Tensor[Tensor-cell2cell]
    Tensor --> TensorExt[Extended Tutorials]
    Compare -->|No| Steady[Steady-state LR Inference]

    %% Multi-modal branch
    Modal -->|Spatial| SMA[Multi-Modal Spatial]
    Modal -->|Non-Spatial| SCMulti[Multi-Modal Single-Cell]

    %% Metabolite sub-branch
    SCMulti --> Metab[Metabolite-mediated CCC]

    %% Links (click events)
    click Inflow "https://liana-py.readthedocs.io/en/latest/notebooks/inflow_score.html"
    click Local "https://liana-py.readthedocs.io/en/latest/notebooks/bivariate.html"
    click Global "https://liana-py.readthedocs.io/en/latest/notebooks/bivariate.html"
    click MISTy "https://liana-py.readthedocs.io/en/latest/notebooks/misty.html"
    click Targeted "https://liana-py.readthedocs.io/en/latest/notebooks/targeted.html"
    click MOFA "https://liana-py.readthedocs.io/en/latest/notebooks/mofatalk.html"
    click Tensor "https://liana-py.readthedocs.io/en/latest/notebooks/liana_c2c.html"
    click TensorExt "https://ccc-protocols.readthedocs.io/en/latest/"
    click Steady "https://liana-py.readthedocs.io/en/latest/notebooks/basic_usage.html"
    click SMA "https://liana-py.readthedocs.io/en/latest/notebooks/sma.html"
    click SCMulti "https://liana-py.readthedocs.io/en/latest/notebooks/sc_multi.html"
    click Metab "https://liana-py.readthedocs.io/en/latest/notebooks/sc_multi.html#metabolite-mediated-ccc-from-transcriptomics-data"
```

## API
For further information please check LIANA's [API documentation](https://liana-py.readthedocs.io/en/latest/api.html).

## Cite LIANA+:

Dimitrov D., Schäfer P.S.L, Farr E., Rodriguez Mier P., Lobentanzer S., Badia-i-Mompel P., Dugourd A., Tanevski J., Ramirez Flores R.O. and Saez-Rodriguez J. LIANA+ provides an all-in-one framework for cell–cell communication inference. Nat Cell Biol (2024). https://doi.org/10.1038/s41556-024-01469-w

Dimitrov, D., Türei, D., Garrido-Rodriguez M., Burmedi P.L., Nagai, J.S., Boys, C., Flores, R.O.R., Kim, H., Szalai, B., Costa, I.G., Valdeolivas, A., Dugourd, A. and Saez-Rodriguez, J. Comparison of methods and resources for cell-cell communication inference from single-cell RNA-Seq data. Nat Commun 13, 3224 (2022). https://doi.org/10.1038/s41467-022-30755-0

Please also consider citing any of the methods and/or resources that were particularly relevant for your research!

[uv]: https://github.com/astral-sh/uv
[scverse discourse]: https://discourse.scverse.org/
[issue tracker]: https://github.com/saezlab/liana-py/issues
[tests]: https://github.com/dbdimitrov/liana-py/actions/workflows/test.yaml
[documentation]: https://liana-py.readthedocs.io
[changelog]: https://liana-py.readthedocs.io/en/latest/release_notes.html
[api documentation]: https://liana-py.readthedocs.io/en/latest/api.html
[pypi]: https://pypi.org/project/liana
