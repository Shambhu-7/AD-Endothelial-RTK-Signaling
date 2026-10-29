# Alzheimer's Disease Brain Endothelial Transcriptomics

## Reproducible receptor tyrosine kinase signaling signatures
## in Alzheimer's disease brain endothelium

Computational analysis of human brain endothelial transcriptomic
data from two independent GEO accessions.

## Datasets

Discovery:
GSE252921
Brain region: Prefrontal cortex
20 donors: 10 Alzheimer's disease, 10 controls

Validation:
GSE327289
Brain region: Inferior temporal gyrus (BA20)
23 donors: 16 Alzheimer's disease, 7 controls

These are separate GEO cohorts. Donor-level independence
across the underlying studies has not been conclusively
established.

## Analysis

- Endothelial pseudobulk aggregation
- Donor-aware differential expression using edgeR
- Age and sex adjustment
- Batch adjustment in the discovery cohort
- Shared-gene cross-cohort comparison
- MSigDB Hallmark and Reactome enrichment
- fgseaMultilevel pathway analysis
- Donor and cell-composition sensitivity analyses

## Main findings

Discovery:
18,321 genes tested; 5 genes at FDR < 0.05.

Validation:
11,029 genes tested; 1 gene at FDR < 0.05.

Shared genes:
10,901.

No gene was FDR-significant in both cohorts.

Pathways tested in both cohorts:
1,168.

One pathway met the same-direction FDR < 0.05
criterion in both cohorts:

REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES

## Reproducibility

The corrected shared-gene pathway analysis was independently
rerun using the released Script 17.

The rerun reproduced all 1,168 discovery and validation
NES and adjusted P values within numerical tolerance,
as well as the identity of the replicated pathway.

## Directory structure

figures/
  Corrected manuscript Figures 2–4

source_data/
  Frozen source data supporting Figures 2–4

scripts/
  R and Python analysis scripts

## Data availability

Original data are available through NCBI GEO:

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE252921

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327289

Raw and processed input matrices are not included
in this GitHub release.

The scripts require the input data to be organized
according to their documented relative paths.

## Interpretation

Transcriptomic enrichment of receptor tyrosine kinase
signaling does not establish RTK protein activation,
blood-brain barrier permeability, or causal mechanisms.

This repository supports computational analysis only.

## Release status

Release candidate v1.0.0.

### Completed verification

- Figure 1: four source-data CSV comparisons passed; regenerated PNG SHA-256 matched.
- Figure 2: three source-data CSV comparisons passed; regenerated PNG SHA-256 matched.
- Figure 3: regenerated PNG SHA-256 matched using frozen pathway inputs.
- Figure 4: regenerated source-data comparison passed; PNG SHA-256 matched.
- Corrected cross-cohort pathway analysis: 1,168 pathway results matched the frozen comparison within the stated numerical tolerance.
- The sole concordantly replicated FDR-significant pathway was REACTOME_SIGNALING_BY_RECEPTOR_TYROSINE_KINASES.

### Reproduction limitations

- The release includes manuscript figures, figure source data, and selected analysis scripts.
- Raw and processed GEO input matrices are not included.
- Some analysis scripts require intermediate files under the original results/ directory.
- Figure-level regeneration has been verified in the original analysis environment.
- Full end-to-end reproduction from downloaded GEO data using only this release has not been demonstrated.
- Independent cohort donor overlap has not been fully resolved.
- GitHub publication, Zenodo archival, and DOI assignment remain pending.
