# pdr-fvm-spp1

Code accompanying a computational reanalysis of proliferative diabetic retinopathy fibrovascular membranes (PDR-FVM), focused on an SPP1-associated cross-state macrophage programme.

## Scope of this repository

This private repository contains only the ten main analysis scripts and documentation. It does **not** contain raw data, intermediate R objects, generated results, figure source files, author information, institutional information, or a software lock file.

The primary public datasets are:

- GSE165784: primary PDR-FVM single-cell dataset.
- GSE102485: exploratory bulk RNA-seq tissue-level comparison.

GSE102485 is affected by disease-status and tissue-source confounding. It must not be interpreted as independent validation or diagnostic validation.

## Current archived-run status

The versions of the ten main scripts archived here were successfully rerun for the analyses supporting the current manuscript. This repository is deliberately a lightweight code archive rather than a self-contained data package.

- Raw public inputs, intermediate R objects, generated results, and figure source files are intentionally excluded from this private archive.
- A new user must obtain the public GEO inputs and create the expected input files before running the workflow independently.
- A future public reproducibility release should add a data-download step, freeze the R and Python environments, and archive the verified runtime record.

## Repository layout

```text
scripts/       Main analysis scripts, executed from the repository root
docs/          Run order, expected outputs, and evidence boundaries
environment/   Historical environment record and future locking instructions
```

## Main analysis order

Run all commands from the repository root. The detailed expected outputs and interpretation boundaries are in [docs/run-order-and-outputs.md](docs/run-order-and-outputs.md).

1. `scripts/01_DR_single_cell_processing_annotation.R`
2. `scripts/02_macrophage_subtype_DEG_GSEA.R`
3. `scripts/03_macrophage_subtype_GO_KEGG.R`
4. `scripts/04_macrophage_subtype_DEG_intersection.R`
5. `scripts/05_bulk_RNA_DESeq2.R`
6. `scripts/06_SPP1_pathway_expression_diagnostic_analysis.R`
7. `scripts/07_SPP1_positive_macrophage_analysis.R`
8. `scripts/08_cell_cell_communication_CellChat.R`
9. `scripts/09_macrophage_pseudotime_analysis.R`
10. `scripts/10_SPP1_transcription_factor_analysis.R`

## Environment

A cleaned historical environment summary is available in [environment/historical-environment.md](environment/historical-environment.md). It is not an executable environment lock and does not prove the software versions used for the earlier result generation. No `renv.lock` is included at this stage.

## Privacy and publication status

The repository intentionally contains no personal names, affiliations, email addresses, local machine paths, author identifiers, or credentials. It is private and has no license. Author, contact, license, citation, release, and DOI metadata must be added only when a public release is approved.
