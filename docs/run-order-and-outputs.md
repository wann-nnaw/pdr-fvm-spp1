# Run order and expected outputs

## Important status note

The mappings below describe the intended relationship between the scripts and the final Figure 1–7 and Supplementary Table S1–S7 materials. They have not yet been confirmed by a clean end-to-end rerun of the copied August 2026 scripts.

| Order | Script | Main purpose | Expected manuscript support |
| --- | --- | --- | --- |
| 01 | `01_DR_single_cell_processing_annotation.R` | Quality control, SCT integration, clustering, primary cell annotation, and the shared annotated object. | Figure 1 and the primary macrophage-state visualisation used in Figure 2; Supplementary Tables S1–S2. |
| 02 | `02_macrophage_subtype_DEG_GSEA.R` | Macrophage-state visualisation, one-versus-rest differential expression, and GSEA. | Figure 2; supporting state-level DEG and enrichment outputs. |
| 03 | `03_macrophage_subtype_GO_KEGG.R` | GO and KEGG over-representation analysis of state-level DEGs. | Figure 3; Supplementary Table S3. |
| 04 | `04_macrophage_subtype_DEG_intersection.R` | Intersection and direction classification across the four macrophage states. | Supplementary Table S4; input for the bulk intersection. |
| 05 | `05_bulk_RNA_DESeq2.R` | Exploratory GSE102485 DESeq2 comparison and overlap with state-level genes. | Supplementary Table S6. |
| 06 | `06_SPP1_pathway_expression_diagnostic_analysis.R` | Exploratory SPP1-related expression and signature analyses. | Supplementary support only; no independent validation or diagnostic claim. |
| 07 | `07_SPP1_positive_macrophage_analysis.R` | SPP1-high versus SPP1-low composition, module scores, DE, and sensitivity analyses. | Figure 5; Supplementary Table S7. |
| 08 | `08_cell_cell_communication_CellChat.R` | CellChat analysis across the primary cell populations. | Figure 4; Supplementary Table S5. Results are predicted cell-cell communication. |
| 09 | `09_macrophage_pseudotime_analysis.R` | Monocle 2 DDRTree pseudotime ordering and sensitivity checks. | Figure 6. This is transcriptional ordering, not lineage or developmental direction. |
| 10 | `10_SPP1_transcription_factor_analysis.R` | Candidate TF, pseudobulk, and pySCENIC regulon analyses. | Figure 7. These are candidate regulatory associations, not direct causal regulation. |

## Required inputs and dependencies

- Step 01 expects GSE165784 input under `data/GSE165784_data/` and writes `rda/UMAP_annotated.rda` for downstream scripts.
- Step 05 expects `data/GSE102485_expressed_gene_reads.txt.gz` and the step 04 intersection output.
- Steps 02, 06–10 depend directly or indirectly on the annotated object created by step 01.
- Step 10 additionally uses the step 09 pseudotime object and pySCENIC outputs. A fresh standard pySCENIC run requires separately obtained cisTarget resources.

## Interpretation boundaries

- SPP1-high is an operational expression-defined group across macrophage states, not a fifth discrete macrophage subtype.
- CellChat results are computational predictions based on transcript and database information.
- Pseudotime results are computational ordering by transcriptional similarity.
- The GSE102485 comparison is exploratory and tissue-level because disease group and tissue source are confounded.
