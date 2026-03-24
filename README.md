# Cell Surface Selectivity Target Explorer

A fast, interactive web application built with FastAPI and vanilla JavaScript to discover, visualize, and analyze cell surface proteins that are highly selectively expressed in specific tumor types relative to normal healthy tissues.

## Overview
Identifying robust cell surface targets is a critical bottleneck in the development of targeted cancer therapeutics (like ADCs or CAR-T therapies). The ideal target is highly expressed in tumors but has minimal or non-existent expression across vital normal organs to prevent on-target, off-tumor toxicity.

This web app integrates over **17,000 RNA-seq samples** from the TCGA (The Cancer Genome Atlas) and GTEx (Genotype-Tissue Expression) databases via the UCSC Xena API. It processes over 60,000 genes to calculate a **Target Selectivity Score** specifically designed to mathematically reward tumor-specific targets.

## Key Features
- **Dynamic Pan-Cancer Boxplots:** View real-time expression distributions of any gene across 32 TCGA primary tumor types, 24 TCGA adjacent normal cohorts, and 51 GTEx completely normal tissues.
- **Selectivity Scoring Metrics:** Ranks over 5,700 Uniprot-verified cell surface proteins based on their targetability.
- **Multi-omics Subtype Splitting:** Dynamically split the TCGA tumor boxes into distinct populations based on:
  - Canonical Molecular Subtypes (e.g., PAM50 in BRCA, IDH status in LGG).
  - Mutational Status (e.g., KRAS, TP53, BRAF, EGFR mutations).
  - Microsatellite Instability (MSI vs MSS) specifically for Colon and Rectal Adenocarcinomas.
- **Sub-second Rendering:** The backend intelligently caches phenotype schemas upon startup and performs surgical HTTP Range requests to the 1.5GB TOIL dataset hosted on Xena via `xenaPython`, eliminating the need to preload entire transcriptome matrices into RAM.

## The Selectivity Score Formula
The backend pipeline evaluated the entire transcriptome line-by-line using the rigorous scoring system below:
1. **Max Tumor Expression ($T_{max}$)**: Maximum median expression across all 32 tumor types.
2. **Max Normal Expression ($N_{max}$)**: Maximum median expression across all 51 GTEx normal tissues.
3. **Tumor vs Normal Delta ($\Delta_{TN}$)**: $\max(0, T_{max} - N_{max})$
4. **Tumor Specificity ($\tau_{tumor}$)**: The tissue-specificity index calculated *only* across the tumor medians. A Tau near 1 signifies expression is highly localized/restricted to just 1 or 2 specific tumor types.
5. **Final Score ($Score$)**: $\Delta_{TN} \times \tau_{tumor}$
   - *Interpretation:* A high Selectivity Score mandates a target must be substantially elevated above its highest normal tissue baseline AND strictly localized to a subset of cancers instead of being uniformly expressed across all of them.
6. **$Score\_T\_AdjN$**: A secondary score isolating strictly the maximum median from the identical *Adjacent Normal* tissue (TCGANormal) instead of the generic GTEx.

## Installation & Local Development

### Prerequisites
- Python 3.8+
- The `surface_selectivity_scores.csv` file (compiled from the raw analysis pipeline) in the root directory.

### Quickstart
1. Clone the repository:
```bash
git clone https://github.com/yourusername/cell-surface-selectivity.git
cd cell-surface-selectivity
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

3. Boot the FastAPI Server:
```bash
uvicorn app:app --reload
```

4. View the App:
Navigate to `http://localhost:8000` in your browser! 
*(Note: The very first time the server spins up, it will take about 20 seconds to completely cache the pan-cancer mutational and phenotype schemas from Xena).*

## Deployment
This application is fully container-ready and configured to be seamlessly deployed via platforms like Render, Railway, or Heroku. The included `Procfile` correctly routes the ASGI server using dynamic port injection. 

Because the architecture streams large arrays dynamically, the web app can run superbly on a standard 512MB RAM Linux container instance.

## Data Sources & Attributions
- Expression Matrices and Phenotype Annotations are programmatically streamed directly from the [UCSC Xena Browser APIs](https://xena.ucsc.edu/).
- **TOIL Dataset:** *TcgaTargetGtex_rsem_gene_tpm*
- **Cell Surface Gene Ontology:** Mined from the Uniprot Subcellular Location database ("Cell membrane/Plasma membrane").
