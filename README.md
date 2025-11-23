# CAR T-cell Atlas for Identifying Patient Response Variability with Single-cell Integration

**Author:** Mohammed Emam Shebl Attia Khattab  
**Supervisor:** Gabriele Sales

---
## Overview
CAR T-cell therapy offers groundbreaking treatment for hematological malignancies like B-cell lymphoma. However, response outcomes vary greatly among patients. This project aims to identify cellular characteristics that distinguish between complete and incomplete responses using an integrated single-cell atlas.

---

## Project Goals
- Explore variability in CAR T-cell therapy outcomes.
- Generate an atlas from **pre- and post-infusion PBMC samples**.
- Use single-cell RNA-seq data to analyze cellular composition and behavior.
- Identify clusters and pathways linked to **patient responses**.

---

## Dataset Summary (all samples included according the project's criteria)
| Time Point        | Sample Count | Cell Count |
|------------------|---------------|------------|
| Pre-infusion     | 24            | ~120,000   |
| Day 7 Post-infusion | 65            | ~260,000   |

---

## Methodology

### 1. Data Collection & Pre-processing
- Pre- and post-infusion PBMC datasets collected from three studies:
#### a. [Post-infusion CAR TReg cells identify patients resistant to CD19-CAR therapy](https://doi.org/10.1038/s41591-022-01960-7) by Good et. al
#### b. [Distinct cellular dynamics associated with response to CAR-T therapy for refractory B cell lymphoma](https://doi.org/10.1038/s41591-022-01959-0) by Haradhvala et. al
#### c. [CAR+ and CAR− T cells share a differentiation trajectory into an NK-like subset after CD19 CAR T cell infusion in patients with B cell malignancies](https://doi.org/10.1038/s41467-023-43656-7) by Louie et. al
- Preprocessed to filter, normalize, and standardize cell-level data.

### 2. Batch Effect Correction
- Applied the Integration tools **Harmony**, **Scanorama**, **CCA**, and **scMerge2**.
- Integration ensures fair comparison across datasets.
![Alt text](/images/preinf_clusters_scanorama)
![Alt text](/images/postinf_clusters_scanorama)


### 3. Cell Type Annotation
- Annotated clusters based on canonical markers.
- Identified key immune populations including CD8+ TEM, CD4+ TCM, NK, MAIT, Mono, etc.

### 4. Composition Analysis
- Applied `sccomp` for statistical testing:
```r
sccomp_result <- object |>
  sccomp_estimate(formula_composition = ~ response + Sex + Age, .sample = id, .cell_group = seurat_clusters) |>
  sccomp_remove_outliers() |>
  sccomp_test()
```
- Selected clusters significantly associated with therapy response.

### 5. Pathway Enrichment
- Conducted GSEA for enriched biological pathways.
- Highlighted immune-related processes in responders.

---

## Key Results

### Pre-infusion
- Enrichment in **leukocyte cell-cell adhesion** and **cytoplasmic translation**.
- Identified cell types (e.g., CD8+ TEM, Proliferating NK, CD14 Mono).

### Day 7 Post-infusion
- Group 1: **OXPHOS**, **Viral Process**, **Adhesion**.
- Group 2: **Regulation of T-cell activation**.
- Identified cell types (e.g., CD8+ naive, CD4+ TCM, dnT, MAIT, ILCs).

---

## Discussion
- Pre-infusion: Cells appear primed, lacking exhaustion phenotypes.
- Post-infusion Group 1: Memory-like metabolic profile for durability.
- Post-infusion Group 2: Balanced activation avoids over-stimulation.

---

## Conclusion
- Persistent responders show unique immune landscapes both pre- and post-infusion.
- Integration-based analysis revealed composition shifts and functional enhancements.
- Suggests that tailored monitoring strategies could predict therapeutic success and inform preemptive measures in cases of likely treatment failure.

---

## Future Directions
- Evaluate viral reactivation markers (e.g., CMV).
- Expand dataset to improve predictive models.
- Examine tumor origin effects on CAR T-cell efficacy.
- Identify biomarkers for therapy optimization.

---
