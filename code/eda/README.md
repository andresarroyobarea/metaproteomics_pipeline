# Metaproteomics EDA Framework

## 1. Overview

This document defines a conceptual framework for exploratory data analysis (EDA) in metaproteomics.

The objective is to systematically characterize data quality, structural properties, abundance patterns, biological signals across multiple biological levels and quantitative metrics.

EDA is organized along three orthogonal analytical dimensions:

1. **Feature-centric analysis** (feature space)
2. **Sample-centric analysis** (sample space)
3. **Condition-centric analysis** (experimental space)

Each dimension is decomposed into conceptual blocks addressing key questions related to quality control (QC), structure, abundance and biology.

This framework is designed to be:

- level-agnostic (peptides, proteins, taxonomy, functions),
- metric-aware (intensity, spectral counts, LFQ-based metrics),
- modular and scalable,
- aligned with preprocessing and filtering logic.

---

## 2. Conceptual Data Types in EDA

Before defining plots, we distinguish three fundamental data objects.

### 2.1 Feature Matrix

Quantitative matrix of features × samples.

- Biological levels:
  - peptide
  - protein
  - taxonomy
  - functional
- Metrics:
  - intensity
  - max_lfq_intensity
  - spectral_count
  - derived metrics (presence/absence, normalized values)

Representations:

- raw values
- normalized values
- presence/absence matrices

---

### 2.2 Feature Metadata

Annotations associated with features.

Examples:

- peptide length, charge state
- protein length, number of unique peptides
- taxonomic ranks
- functional categories

---

### 2.3 Sample Metadata

Annotations associated with samples.

Examples:

- condition
- batch
- clinical variables
- inclusion/exclusion flags

---

## 3. Feature-centric Analysis

**Central question:**  
👉 How do features behave in terms of abundance, prevalence, stability and biological properties?

This dimension describes the structure of the feature space independently of experimental conditions.

---

### 3.1 Quality Control (QC) at Feature Level

#### 3.1.1 Feature Abundance Distribution

**Objective:**  
Evaluate the global distribution of feature abundance.

**Plots:**

- histogram / density plot of metrics
- log-scale distributions

**Example:**

- X: log10(intensity or spectral count)
- Y: frequency of features

---

#### 3.1.2 Feature Dominance / Heterogeneity

**Objective:**  
Assess whether the signal is dominated by a small subset of features.

**Plots:**

- rank-abundance curve (log scale)
- cumulative abundance curve
- Lorenz curve (Gini-like)

**Example:**

- X: feature rank (sorted by abundance)
- Y: cumulative abundance (% of total signal)

---

### 3.2 Prevalence Structure

#### 3.2.1 Feature Prevalence Distribution

**Objective:**  
Assess how frequently features are detected across samples.

**Plots:**

- histogram or barplot of prevalence

**Example:**

- X: number of samples where a feature is present
- Y: number or percentage of features

---

#### 3.2.2 Prevalence vs Abundance

**Objective:**  
Evaluate the relationship between abundance and prevalence.

**Plots:**

- scatter plot / hexbin plot

**Example:**

- X: mean intensity
- Y: prevalence (% samples)
- color: metric type (intensity vs spectral count)

---

#### 3.2.3 Feature Consistency

**Objective:**  
Assess variability of features across samples.

**Plots:**

- histogram of coefficient of variation (CV)
- mean–variance plot

**Example:**

- X: mean intensity
- Y: CV
- color: metric type

---

### 3.3 Filtering Impact Analysis

**Objective:**  
Quantify the impact of preprocessing filters on the feature space.

**Plots:**

- before vs after filtering barplots
- stacked barplots
- flowcharts / Sankey-like plots

**Examples of filters:**

- non-human filter
- unique peptide/protein filter
- minimum prevalence filter
- condition-specific filters

**Output:**

- number of features retained at each filtering step

---

### 3.4 Biology-specific Feature Properties

#### Peptide-level

- peptide length distribution
- charge state distribution

#### Protein-level

- protein length distribution
- number of unique peptides per protein

#### Taxonomy-level

- taxonomic rank distribution
- relative abundance by rank

---

## 4. Sample-centric Analysis

**Central question:**  
How do samples behave in terms of signal, complexity, quality and structural similarity?

This dimension characterizes sample-level properties independently of biological interpretation.

---

### 4.1 Quality Control (QC) at Sample Level

#### 4.1.1 Total Abundance per Sample

**Plot:** barplot

**Example:**

- X: sample ID
- Y: total intensity / total spectral count

---

#### 4.1.2 Feature Counts per Sample

**Plot:** barplot

**Example:**

- number of detected features per sample
- comparison across metrics

---

#### 4.1.3 Missingness per Sample

**Plots:**

- barplot of missing values
- missingness heatmap

**Example:**

- % of zero or NA values per sample

---

#### 4.1.4 Sample Complexity

**Objective:**  
Assess whether samples have comparable complexity.

**Plots:**

- rarefaction curves
- cumulative intensity curves

**Example:**

- X: number of features
- Y: cumulative intensity (%)

---

#### 4.1.5 Sample Dominance

**Objective:**  
Assess whether a sample is dominated by few features.

**Plots:**

- top-N contribution plot
- cumulative abundance curves per sample

**Example:**

- % of total intensity explained by top 10 / top 50 features

---

### 4.2 Abundance Distributions Across Samples

**Plots:**

- boxplots / violin plots

**Example:**

- distribution of intensities per sample

---

### 4.3 Sample Structure

#### 4.3.1 Similarity Patterns

**Plots:**

- feature × sample heatmaps
- distance heatmaps

---

#### 4.3.2 Ordination / Dimensionality Reduction

**Plots:**

- PCA
- PCoA / NMDS
- t-SNE / UMAP

**Note:**  
Ordination should be performed on a subset of features (e.g. filtered or high-variance features).

---

#### 4.3.3 Clustering

**Plots:**

- hierarchical clustering dendrograms

---

#### 4.3.4 Sample Correlation

**Objective:**  
Evaluate coherence of technical and biological replicates.

**Plots:**

- sample–sample correlation heatmaps
- pairwise scatterplots

---

### 4.4 Normalization Diagnostics

**Objective:**  
Assess the impact of normalization on data structure.

**Plots:**

- PCA pre vs post normalization
- density plots pre vs post normalization
- boxplots pre vs post normalization

---

## 5. Condition-centric Analysis

**Central question:**  
👉 How do experimental conditions differ and what biological patterns emerge?

This dimension focuses on experimental design and biological contrasts.

---

### 5.1 QC and Bias Detection by Condition

#### 5.1.1 Metric Distribution by Condition

**Plots:**

- boxplots
- density plots

---

#### 5.1.2 Sample Balance

**Plots:**

- barplot of sample counts per condition

---

### 5.2 Abundance Patterns by Condition

#### 5.2.1 Metric Abundance

**Plot:** boxplot by condition

---

#### 5.2.2 Feature Counts by Condition

**Plot:** barplot

---

#### 5.2.3 Condition Effect Size

**Plots:**

- fold-change distributions
- log2FC density plots

*(typically more relevant after differential analysis)*

---

### 5.3 Feature Overlap Between Conditions

#### 5.3.1 Shared vs Specific Features

**Question:**  
Are there condition-specific features?

**Plots:**

- Venn diagrams (few conditions)
- UpSet plots (recommended)

---

### 5.4 Condition Structure

#### 5.4.1 Ordination by Condition

**Plots:**

- PCA / PCoA colored by condition

---

#### 5.4.2 Drivers of Separation

**Objective:**  
Identify features driving condition differences.

**Plots:**

- PCA loadings
- top discriminant features plots

---

#### 5.4.3 Beta Diversity

**Plots:**

- distance matrices
- ordination plots based on distance metrics

---

## 6. Conceptual Notes on Metrics in Metaproteomics

In this framework, abundance is treated as a multi-metric concept:

- intensity (MS1-based abundance)
- max LFQ intensity (normalized intensity-based abundance)
- spectral count (MS2-based detection/abundance proxy)

Presence/absence is defined as a metric-dependent derived layer.

Not all metrics are available at all biological levels, and this is explicitly acknowledged in the EDA design.

---

## 7. Design Philosophy

This EDA framework is designed to be:

- level-agnostic (peptides, proteins, taxonomy, functions)
- metric-aware (intensity, spectral counts, LFQ)
- modular (reusable plotting functions)
- scalable (from exploratory plots to publication-ready figures)
- pipeline-compatible (aligned with preprocessing and filtering logic)
---
