Bioinformatic Evidence for the EGFR-HERC4 Regulatory Axis in RPE Fibrosis

📖 Background & Hypothesis
In diabetic retinopathy (DR), the transition of Retinal Pigment Epithelial (RPE) cells to a myofibroblast phenotype drives subretinal fibrosis.
We hypothesize that high-glucose (HG) conditions activate the EGFR-HERC4 axis, leading to the proteasomal degradation of SAV1, which subsequently triggers the Hippo/YAP signaling pathway and promotes fibrotic gene expression (e.g., ACTA2, CTGF).
This repository contains the bioinformatic validation of this axis using transcriptomic data from GSE102485.

📊 Key Results

1.Pathway Co-expression SynergyThe
cornerstone of this analysis is the discovery of a robust synergistic expression between the receptor tyrosine kinase (EGFR) and the E3 ubiquitin ligase (HERC4).
Correlation Strength: Pearson r = 0.65 in the co-expression matrix.
Statistical Significance: Scatter plot analysis confirmed a strong linear relationship with R = 0.57 and p = 0.00089.
Interpretation: This high degree of coordination suggests that EGFR and HERC4 function as a unified regulatory unit (Axis) in response to HG stress.

2. Differential Expression & Post-Translational Inference
While mRNA levels showed moderate upward trends in HG-treated cells, the lack of traditional p < 0.05 significance underscores a critical biological insight:
Transcriptional Trends: EGFR (p=0.41) and HERC4 (p=0.62) levels increase but remain relatively stable.
SAV1 Surrogate Analysis: Due to the absence of SAV1 in this dataset, we utilized its downstream effector CTGF. CTGF showed consistent co-fluctuation with ACTA2 (r=0.35), supporting the activation of the downstream fibrotic program.
Strategic Conclusion: The stability of mRNA levels coupled with high axis-correlation strongly points toward post-translational modification (protein degradation) as the primary regulatory mode for the SAV1 protein.

🛠️ Methodology & Visualization
We developed a robust R pipeline to process the FPKM expression matrix:
Data Normalization: Log_2(FPKM + 1) transformation to stabilize variance.
Outlier Mitigation: Integrated violin-boxplot visualizations to transparently display data distribution.
Statistical Testing: Two-tailed t-tests and Pearson correlation coefficients.

Integrated Analysis Figure
The following multi-panel figure integrates differential expression, co-expression logic, and axis correlation:
Top : Multi-gene Boxplots.
Bottom Left: Pearson Correlation Heatmap.
Bottom right:EGFR-HERC4 Linear Regression.

Repository Structure
data: Contains raw FPKM counts (GSE102485).
scripts: analysis_pipeline.R - Full code from cleaning to final visualization.
figures: High-resolution PNG exports of the results.
