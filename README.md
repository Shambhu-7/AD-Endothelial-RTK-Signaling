# Alzheimer's Disease Neurovascular Biology: Paper 4

## Cross-cohort analysis of human brain endothelial transcriptomics

This repository contains computational analyses investigating endothelial transcriptional changes associated with Alzheimer's disease (AD), with particular attention to receptor tyrosine kinase (RTK) signaling, pathway reproducibility, and neurovascular research priorities.

## Datasets

| Role | GEO accession | Donors | AD | Controls |
|---|---|---:|---:|---:|
| Discovery | GSE252921 | 20 | 10 | 10 |
| Validation | GSE327289 | 23 | 16 | 7 |

Donor-level independence across the underlying studies has not been conclusively established. The cohort counts should not be interpreted as proof of independent sampling.

## Analysis workflow

- Endothelial pseudobulk expression analysis
- Donor-aware differential expression using edgeR
- Covariate adjustment and cell-composition sensitivity analyses
- Cross-cohort gene-effect comparison
- Hallmark and Reactome pathway enrichment
- Pathway replication and robustness assessment
- Exploratory computational network and machine-learning analyses

## Principal findings

- Discovery analysis: 5 genes passed FDR < 0.05 among 18,321 tested.
- Validation analysis: 1 gene passed FDR < 0.05 among 11,029 tested.
- No shared gene was FDR-significant in both cohorts.
- Of 1,168 compared pathways, RTK signaling was the only pathway meeting the same-direction FDR < 0.05 replication criterion in both cohorts.

These findings support a reproducible RTK-signaling enrichment association. They do not establish RTK protein activation, blood–brain barrier dysfunction, or causal disease mechanisms.

## Repository structure

- `release_v1.0.0/`: versioned figures, source data, selected analysis scripts, environment records, and release documentation.
- `scripts/`: working analysis scripts, including exploratory analyses.
- `results/`: intermediate analysis outputs and audit results; not all outputs are part of the public release.
- `data/`: local data files, excluded from Git tracking.

## Reproducing the released figures

Use Python 3 with the packages listed in `requirements.txt`.

Activate the project's environment and run:

```bash
python -m pip install -r requirements.txt
python release_v1.0.0/scripts/regenerate_figures_from_release.py
```

The script reconstructs diagnostic figures from the released source-data CSVs. These reconstructions are not pixel-identical replacements for the original manuscript figures.

The complete original analysis workflow has not yet been demonstrated from raw GEO downloads using only the public release. Some analysis scripts require intermediate files from the original working environment.

## Data availability

Original datasets are available from NCBI Gene Expression Omnibus:

- [GSE252921](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252921)
- [GSE327289](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327289)

Raw and processed GEO input matrices are not included in the versioned release.

## Current limitations

- Donor-level independence between cohorts remains unresolved.
- Full end-to-end reproducibility from raw public data has not been demonstrated.
- Exploratory AI and network analyses require independent evaluation before mechanistic conclusions can be drawn.
- Zenodo archival and DOI assignment remain pending.

## Scope

This repository supports computational research into AD-associated endothelial transcriptional changes. Experimental validation is required to establish functional effects on the neurovascular unit and blood–brain barrier.

