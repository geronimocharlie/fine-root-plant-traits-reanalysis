

### **1 Imputation uncertainty not propagated**

* The study imputes missing trait values (particularly for fine-root traits) using a single run of the *missForest* algorithm with phylogenetic predictors.
* However, uncertainty from this imputation is **not quantified or propagated** into downstream analyses (PCA, trait probability densities, redundancy metrics, group comparisons).
* This omission may underestimate variance and overstate precision.
* **Suggestion:**

  * Conduct **multiple imputations** (same method, different random seeds) to assess stability of PCA loadings and group dissimilarities; pool results (Rubin’s rules).
  * Alternatively, use a **parametric bootstrap** drawing from imputation error distributions (out-of-bag MSE from *missForest*) to propagate uncertainty through the multivariate analyses.

---

### **2 Heterogeneity in belowground data (field vs pot studies)**

* The fine-root trait dataset (RES) integrates data from studies conducted under **diverse experimental conditions** (field, pot, hydroponic).
* While the authors show that mean trait values correlate strongly (>0.97) between field-only and full data, this does **not ensure covariance-level consistency** — the multivariate relationships that define the PCA could still differ.
* **Why this matters:** PCA depends on how traits *co-vary*, not just on their means; measurement context can alter those relationships.
* **Suggestions:**

  * Re-run the **entire PCA and TPD analyses** on **field-only data** and compare to the full dataset using a **Procrustes or RV test** to quantify stability of the trait space.
  * As an intermediate step, test a **weighted PCA**, where each species is weighted by data reliability (e.g., proportion of field measurements, number of independent studies).
  * In the kernel density (TPD) analyses, inflate bandwidths (smoother densities) for species derived mainly from pot studies to down-weight uncertain data.

---

### **3 Inconsistent definition of “fine roots”**

* Studies contributing to the RES dataset use **different definitions** of fine roots — either **diameter-based** (e.g. roots ≤ 2 mm) or **order-based** (e.g. first 1–3 branching orders).
* These criteria capture **different functional segments** of the root system: diameter-based samples can include transport roots, while order-based samples focus on absorptive roots.
* This inconsistency means that “fine-root” traits may not refer to the same biological structures across species or studies, creating functional noise.
* **Suggestions:**

  * Quantify how many species are represented by **order-based vs diameter-based definitions** (and those with mixed sources).
  * Re-run the PCA and trait-space analyses after **excluding mixed-definition species**, or separately for each definition type.
  * If results diverge, this should be discussed explicitly as a methodological limitation of global root trait synthesis.

---

### **4 Potential artifact from varimax rotation**

* The use of **varimax rotation** in PCA simplifies interpretation by maximizing simple structure, but it can **mathematically enforce trait separation** between domains.
* Extended Data Table 2 shows that before rotation, PC1 blended **size and root traits**, whereas after rotation, above- and belowground traits fall on distinct, nearly orthogonal axes — a non-trivial conceptual change.
* Although the Procrustes correlation between rotated and unrotated spaces is high (r ≈ 0.99), that metric hides a **qualitative re-partitioning of variance** that affects biological interpretation.
* **Suggestions:**

  * Present both the **rotated and unrotated** PCA results in the main text or supplement, clearly explaining differences in trait loadings and interpretation.
  * Discuss that rotation improves readability but may **overstate independence** of root and leaf domains.
  * Optionally, complement PCA with an **ordination method that doesn’t assume linear orthogonality** (e.g., NMDS or factor analysis without rotation) to verify structural robustness.

