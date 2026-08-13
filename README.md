# endometriosis_study
Exploration of the single-cell atlas data from Marečková et al.

## Conda environment

Create the project environment from the repository root:

```bash
conda env create --file environment.yml
```

Register the environment as a Jupyter kernel:

```bash
conda run --name endometriosis-study \
  python -m ipykernel install --user \
  --name endometriosis-study \
  --display-name "Endometriosis study"
```

In JupyterLab or an IDE, open
`notebooks/endometriosis_scRNAseq_analysis.ipynb` and select the
**Endometriosis study** kernel.
