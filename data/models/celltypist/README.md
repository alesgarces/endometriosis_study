# Human Endometrium Atlas CellTypist model

The official `Human_Endometrium_Atlas.pkl` model was downloaded through CellTypist 1.7.1 on 2026-08-16. CellTypist describes it as 36 endometrial cell types integrated from seven datasets across the menstrual cycle and cites [Marečková *et al.* (2024)](https://doi.org/10.1038/s41588-024-01873-w) as its source.

| File | SHA-256 |
| --- | --- |
| `data/models/Human_Endometrium_Atlas.pkl` | `daafe27a9c444fd2939b4d296dbc24e61223336b85278f9ae0dc95d9af500ede` |

The nested `data/models` directories follow CellTypist's standard cache layout under the project-local `CELLTYPIST_FOLDER`. The notebook downloads the model only when it is absent.

Because this project's cells contributed to the source atlas, predictions reproduce the paper's annotation framework but do not constitute independent validation. Labels still require marker, cluster, donor, library, condition, and doublet checks before biological interpretation.
