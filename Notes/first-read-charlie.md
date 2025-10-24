### **1. Purpose of the paper**

The paper aims to **link aboveground and belowground plant traits** to understand whether plants form a *coordinated “whole-plant economics spectrum”* — meaning: do fast-growing (resource-acquisitive) leaves go with fast, acquisitive roots?

Past studies mostly focused on leaves and stems, but roots were missing or poorly integrated. This study fills that gap.

---

### **2. Data used**

The study merges two major *global plant trait databases*:

- **GSPFF** (Global Spectrum of Plant Form and Function):
    
    A dataset of **six key aboveground traits** measured across ~46,000 species.
    
    Traits include things like **plant height, seed mass, leaf area, leaf thickness, leaf nitrogen, and leaf mass per area (LMA)** — i.e., the ones describing size and “leaf economics.”
    
- **RES** (Fine-Root Economics Space):
    
    A newer dataset of **four fine-root traits** from ~1,500 species.
    
    Traits include **specific root length (SRL), root diameter (D), root tissue density (RTD), and root nitrogen concentration (RNC)** — describing root construction and nutrient strategies.
    

After harmonizing taxonomy and imputing missing data, the authors obtained **301 species with complete measurements** of all 10 traits and **1,218 species with imputed data**.

---

### **3. What combining them means**

By combining GSPFF and RES traits, they could create a **joint “whole-plant trait space”** — a multivariate map of plant strategies spanning leaves *and* roots.

Then they ask:

- How many independent axes (gradients) describe global plant strategies?
- Are leaf and root economics tightly linked or largely independent?
- Do growth forms (woody/herbaceous) or biomes occupy different parts of that space?

---

### **4. Main finding (in simple terms)**

The combined data show **four main axes**:

- Two aboveground ones: *plant size* and *leaf economics*
- Two belowground ones: *root tissue economics* and *root–fungal interaction strategies*

These sets are **partly coordinated but not tightly aligned** — meaning above- and belowground strategies evolve somewhat independently.

---

Would you like me to make a quick diagram (text or visual) showing how GSPFF and RES traits align in those four axes? It can really help make the structure intuitive.

# Summary of contribution

The paper integrates six canonical aboveground traits (GSPFF) with four fine-root traits (RES) to build a joint trait space. Using PCA with varimax rotation on 301 species with complete data—and an imputed set of 1,218 species—the authors show that four orthogonal axes are required: two “aboveground” axes (size and leaf economics) and two “fine-root” axes (root tissue conservation and plant–fungal interaction gradient). The aboveground vs fine-root planes display distinct aggregation patterns (two hotspots aboveground; one hotspot belowground) and markedly different differentiation across growth forms, families, and biomes. The analysis suggests partial whole-plant coordination (correlated dissimilarities among groups) but limited one-to-one alignment of leaf and root economics.

# Major strengths

- **Clear dimensional result with careful cross-checks.** The “four-axis” conclusion is well supported by (i) PCA on complete data; (ii) varimax rotation to improve interpretability; and (iii) Procrustes tests linking the joint space to separate aboveground and root spaces (r ≈ 0.99 and 0.98, respectively), while showing weak direct correspondence between leaf and root spaces.
- **Thoughtful null‐modeling of trait space occupancy.** The TPD/kde approach plus multivariate-normal null profiles gives interpretable evidence for clumping and lower divergence in the fine-root plane vs aboveground.
- **Group-level ecology made quantitative.** Woodiness explains ~36.5% of aboveground variation but ~0.4% on the fine-root plane; families/biomes separate much more aboveground; and redundancy–productivity associations diverge above vs belowground. These are biologically meaningful and statistically explicit patterns.
- **Transparent data processing and availability.** The data sources (TRY, GRooT), taxonomy harmonization, GBIF/biome assignment, and code/data availability via Figshare are documented, enabling reproducibility.

# Points needing clarification or additional analysis

1. **Trait imputation uncertainty is not propagated.**
    
    The missForest + phylogenetic eigenvector imputation is reasonable, but subsequent analyses treat imputed values as known. Please quantify how imputation affects (a) axis loadings, (b) group dissimilarities, and (c) redundancy metrics. Two concrete options:
    
    - Repeat all key statistics on **multiple imputed datasets** (m ≈ 20) and pool estimates (Rubin’s rules); or
    - **Parametric bootstrap**: draw from imputation error distributions (out-of-bag MSE from missForest) to propagate uncertainty into PCA, TPD, PERMANOVA, and Mantel tests.
2. **Heterogeneity in measurement conditions and designs.**
    
    You model study design and publication with mixed models for each root trait and show high correlations field-only vs field+pot (>0.97). That’s helpful, but key multivariate results (PCA/rotation; TPD surfaces) still inherit residual study-level error. Consider:
    
    - **Sensitivity**: re-run the joint PCA/varimax and TPD analyses on field-only roots (you provide trait-wise robustness; extend to the **full pipeline**).
    - **Study-weighted kernels**: inflate bandwidths (or down-weight densities) for species whose trait values come from few studies/sites, to temper overconfidence where sampling is sparse.
3. **Root definition and within-species variance.**
    
    Fine roots are compiled across definitions (≤2 mm vs orders 1–3). Even though you prefer minimum order / diameter, species means may still blend distinct root orders and mycorrhizal zones. Please:
    
    - Report how many species mix **order-based** and **diameter-based** records.
    - Add a sensitivity where you **exclude mixed-definition species** (or use only order-defined records) and show the effect on the C3–C4 plane and on the “single hotspot” result.
4. **Biome assignment via GBIF occurrences can bias climate placement.**
    
    Occurrence density varies massively (1 to 906,097 records per species). This can over-represent well-sampled regions. Suggested checks:
    
    - **Equal-weight species** when averaging climate per biome (i.e., average species-level climate means, not raw occurrence means).
    - **Spatial thinning** of occurrences before climate extraction.
    - Show robustness of biome dissimilarities to these alternatives.
5. **NPP estimation and redundancy relationships.**
    
    Redundancy aboveground ↑ with aboveground NPP, while fine-root redundancy ↓ with belowground NPP (estimated from precipitation). Because NPP is modeled from MAP, this may conflate precipitation with productivity. Consider repeating with a **gridded NPP product** (e.g., MODIS or FLUX-informed) or, at minimum, add a sensitivity using **temperature + precipitation** together to compute NPP surrogates. Also report **partial regressions** (redundancy ~ NPP controlling for biome identity) to guard against small-n leverage (n=8 biomes).
    
6. **Ordination/rotation choices.**
    
    PCA on log-scaled traits is fine for global structure, but species packing and non-Gaussian trait distributions can warp Euclidean distances. Helpful additions:
    
    - Confirm axis structure with **robust PCA** (e.g., projection-pursuit) and/or **nonlinear ordination** (NMDS/GLLVM latent variables) for the 301-species set.
    - Provide **bootstrap axis stability** (loadings and species scores) to quantify uncertainty in Fig. 1 spaces.
7. **Phylogenetic signal and group effects.**
    
    You use phylogeny during imputation and run PERMANOVA for group separations. It would strengthen interpretation to quantify **phylogenetic signal** (Pagel’s λ or Blomberg’s K) on each axis, and to report **phylogenetically informed PERMANOVA** or distance-based PGLS to test whether family/woody vs herbaceous differences persist after controlling for shared ancestry.
    
8. **Kernel density and bandwidth choices.**
    
    The TPD approach uses plug-in bandwidths (Hpi). In high-redundancy settings, bandwidth strongly shapes hotspot topology. Please add:
    
    - A sensitivity panel with **(i)** Hpi, **(ii)** cross-validated bandwidth, **(iii)** a fixed scalar multiple of the pooled covariance, showing that the “one hotspot belowground” conclusion is invariant.
9. **Scale and missing belowground size traits.**
    
    You note the absence of root size/allocation traits as a likely reason for weaker structuring belowground. This is an important limitation. If any partial data exist (rooting depth, coarse-root biomass, clonal spread), a **reduced-trait analysis** could at least test whether adding any belowground size proxy strengthens separation—else prominently flag this as a priority for future data mobilization.
    

# Minor comments & presentation

- **Figure clarity:** In Fig. 1, annotate axes with dominant traits (“Size/leaf economics” and “Root tissue conservation/plant–fungal interactions”) in the panel titles to aid readers skimming the figures. The current arrow legend is good but small.
- **Reporting PERMANOVA:** Provide **R² and PERMDISP** checks to ensure that differences aren’t driven by dispersion heterogeneity among groups.
- **Thresholding choice:** You consistently use the 0.99 TPD quantile; add a supplementary figure with **0.95** to show robustness of overlap/dissimilarity patterns.
- **Units and transformations:** Early in Methods, list **all trait units post-transform** (log10, centered/scaled) once, so readers can quickly interpret loadings across datasets.
- **Terminology:** Where you state “virtually indistinguishable” root syndromes for woody vs herbaceous species, include the **effect size and CI** (e.g., R² ≈ 0.4%) inline to quantify “virtually.”

# 

# Details on each point of critique

## 1 Imputation

---

### **What “Trait imputation uncertainty” means here**

Yes — the authors had **missing trait values** (especially for root traits), and to avoid losing most species, they **imputed** them using a statistical model.

They used **missForest** (a random-forest-based imputation algorithm), augmented with **phylogenetic eigenvectors** as predictors.

This is a sensible method, but — crucially — they **ran it only once**, producing a **single “best guess” dataset**.

They then treated the imputed values as if they were real, known measurements — ignoring the uncertainty that comes from imputation.

---

### **Why this matters**

Imputation uncertainty can influence:

- PCA loadings and axis structure,
- Group differences (e.g., woody vs herbaceous),
- Metrics like trait space redundancy or dissimilarity.

If the imputation introduces bias or overconfidence, it could make those multivariate results look more precise than they really are.

---

### **What they *should* do (two main approaches)**

### **A. Multiple imputation (preferred and standard)**

- Run **the same imputation method** (e.g. missForest with phylogenetic predictors) **multiple times** — say, 20 times — with slightly different random seeds.
- Each run produces one *plausible* completed dataset.
- Then, re-run the analyses (PCA, TPD, etc.) on each dataset separately, and **pool the estimates** (means, variances, or significance) using **Rubin’s rules** or simple averaging.

This quantifies how sensitive the results are to imputation randomness.

They could also test *different imputation methods* (e.g., Bayesian PCA, random forest, kNN) as a **robustness check**, but the key point is *multiple runs*, not just one.

---

### **B. Parametric bootstrap approach**

This approach keeps the same imputation model (missForest), but uses its **out-of-bag prediction errors** to simulate variability:

- For each imputed value, estimate its **prediction error variance** from missForest.
- Generate **bootstrap datasets** by adding random noise (drawn from those error distributions) to the imputed values.
- Re-run the PCA and other analyses many times to assess how much the conclusions vary due to imputation uncertainty.

So yes — in this bootstrap scenario, you **stick to the same imputation method**, but **sample from its uncertainty** rather than re-fitting new imputations from scratch.

---

### **In short**

> They imputed missing traits once with missForest, didn’t propagate that uncertainty, and treated imputed values as fixed.
> 
> 
> A stronger approach would be to run **multiple imputations** (same method, different random seeds) or to **bootstrap imputation error** to see how stable their results really are.
> 

---

## 2 Heterogenity

---

### **1. Where the heterogeneity comes from**

Yes — the heterogeneity is **only (or mostly) in the root dataset (RES)**, because it’s compiled from many independent studies.

Each study may differ in:

- **Growth environment** (field vs pot),
- **Measurement protocol** (definition of root order, washing method, drying time, etc.),
- **Sampling depth or root type.**

This means that, even after unit harmonization, some residual variation reflects *methodological differences*, not biology.

---

### **2. What the authors did**

They modeled each root trait separately using **mixed models** like:

[

\text{Trait} = \text{fixed effects} + (1 | \text{StudyID}) + (1 | \text{Design})

]

and found that **species’ mean trait values** from the *field-only subset* and the *full dataset* were very highly correlated (>0.97).

That’s reassuring for single traits — it means the *average* effect of study design doesn’t strongly bias the species mean for each trait.

---

### **3. Why this doesn’t guarantee the PCA (trait space) is unaffected**

Good question — here’s the subtle but important part:

Even if each trait alone is almost unchanged, **the way traits *covary* across species** can still differ.

PCA is based on **the covariance (or correlation) matrix** of all traits.

So, the relative structure of correlations matters more than individual trait means.

Example:

- Suppose trait A and trait B both shift slightly due to measurement design, but not equally.
- The *direction* of that shift may change the correlation between A and B — and hence rotate the PCA axes.

So, you can have **high pairwise trait correlations between datasets (>0.97)** but still **meaningfully different multivariate structure**, especially in higher dimensions.

That’s why just comparing field-only vs full data for single traits isn’t enough — the *interaction structure* of the traits (the PCA loadings and scores) could still change.

---

### **4. What to do about it**

Yes — your proposed fix is exactly right:

> Run the PCA (and the subsequent trait-space analyses) twice:
> 
> 1. Using the full root dataset (field + pot),
> 2. Using only field-measured root traits.

Then compare the resulting **axes (loadings)** and **species positions (scores)** using a **Procrustes test or RV coefficient**.

If the alignment is very high (r > 0.95), the conclusion is robust.

If not, it shows that experimental heterogeneity affects the derived trait space.

---

### **5. Optional improvement: weighting by data reliability**

This is a more nuanced idea — here are two ways it could be implemented:

### **A. Weighted PCA**

Each species’ contribution to the PCA could be **weighted** by a measure of data reliability, e.g.:

- Number of independent studies that measured that species’ traits,
- Fraction of traits directly measured vs imputed,
- Whether the trait values came from field or experimental conditions.

Technically, this can be done by:

- Using a **weighted covariance matrix**,
    
    (\Sigma_w = (X^T W X) / \sum w_i),
    
    where (W) is a diagonal matrix of species weights.
    
- Or by **duplicating** high-quality species proportionally in the dataset (a rough but simple approach).

### **B. Adjust kernel density estimation (TPD)**

When estimating the **trait probability distributions (TPD)** for species groups, they could:

- Increase the **bandwidth** (smoothing) for species with low data reliability (so their distributions are more diffuse), and
- Decrease it (sharpen) for well-sampled species.

That way, uncertain or heterogeneous data don’t artificially create sharp “hotspots” in the trait space.

---

### **6. Summary interpretation**

| Step | Concept | Explanation |
| --- | --- | --- |
| **Heterogeneity** | Root data come from different experimental conditions. | Measurement differences can distort trait relationships. |
| **What they did** | Checked single-trait robustness (field vs full). | High correlations (>0.97) show stable averages. |
| **What they missed** | Covariance-level effects (trait interactions). | PCA could still differ because relationships among traits may shift. |
| **Fix 1 (essential)** | Re-run PCA on field-only vs full datasets and compare via Procrustes. | Tests if overall trait space structure is stable. |
| **Fix 2 (optional)** | Weighted PCA or kernel weighting. | Down-weight or smooth uncertain data to reduce bias. |

---

## 3 On defining Roots

---

### **1. Two common ways to define “fine roots”**

When researchers measure root traits, they need to decide *which roots* to sample — because a single plant has roots of very different thicknesses and functions.

There are **two main conventions** used across studies:

### **A. Diameter-based definition**

- Fine roots are defined as **roots thinner than a threshold**, typically **≤ 2 mm** in diameter.
- It’s simple to apply, but the problem is that the **functional meaning** of “2 mm” varies among species.
    - In woody plants, roots < 2 mm can include both absorptive and transport roots.
    - In herbs, almost all roots might be < 2 mm.

### **B. Order-based definition**

- Fine roots are defined by their **position (branching order)** in the root system, starting from the tips.
    - **1st-order roots** are the youngest, most distal roots that absorb water and nutrients.
    - Higher orders (2nd, 3rd, etc.) are thicker and more transport-oriented.
- Researchers often sample **orders 1–3** as “fine roots.”

---

### **2. Why this matters**

These two approaches are **not equivalent**, and they can produce quite different average trait values:

- Diameter-based sampling can include thicker, more structural roots.
- Order-based sampling targets only the active, absorptive roots with high specific root length and nitrogen.

If some studies used the diameter criterion and others used branching order, the “fine-root traits” may not refer to exactly the same biological entities.

---

### **3. What happens in this paper**

In the **RES dataset**, the authors combined root data from many studies.

- Some studies reported **order-based** measurements, others **diameter-based**, and for some, it’s not even clear.
- When harmonizing the data, they likely chose (where possible) the **smallest available order or diameter** to represent each species’ fine roots.

This mixing means that, for some species, the recorded root traits might come from **different types of roots** depending on the study.

So yes — in a way, some species could fall into both “fine” and “less fine” categories depending on which data source was used.

---

### **4. Why it’s a problem**

This inconsistency introduces **methodological heterogeneity** again — but specifically in *what is being measured*.

It could blur the differences between species or distort trait correlations, because some “fine roots” are truly absorptive while others are partly structural.

---

### **5. What the critique suggests**

To check robustness:

- Identify species for which data come from **mixed definitions** (both order-based and diameter-based studies).
- Run the main analyses (PCA, trait space mapping) **excluding** those mixed-definition species, or **using only order-defined records**.
- See whether the main structure (the four trait axes, the fine-root hotspot) holds.

If it’s stable, great — if not, it shows that definitional inconsistency affects the inferred “root economic spectrum.”

---

### **In short**

> “Order-based” and “diameter-based” refer to two different ways of defining fine roots.
> 
> 
> Mixing them means some species’ root traits might not represent the same biological level of the root system.
> 
> The critique asks for a sensitivity test excluding or separating these definitions, to ensure that the differences in root trait structure aren’t an artifact of inconsistent fine-root definitions.
> 

---

## 4 GBIF biome representation

---

### **1. Where GBIF fits into the data**

GBIF (the **Global Biodiversity Information Facility**) provides **species occurrence records** — basically, where on Earth each species has been observed.

In this study, GBIF is **not** the source of the trait data themselves.

Instead, it is used to:

- **Map the distribution of each species** (using latitude/longitude records),
- **Extract associated climate variables** (e.g. temperature, precipitation) from those locations, and
- **Assign each species to one or more biomes** (e.g. tropical forest, temperate grassland, boreal, desert, etc.) based on those coordinates.

So GBIF provides the **spatial and environmental context** for the species that have trait data.

---

### **2. Why this can be a problem (sampling bias)**

GBIF data are **unevenly sampled across the globe**.

Some regions (e.g. Europe, North America) have millions of occurrence records per species; others (e.g. tropical Africa, parts of Asia) have very few.

That means:

- For species that occur globally, **the climate/biome mean** might be skewed toward the best-sampled regions.
- For species with only a few records, a single cluster of points could misrepresent their true climate range.

This bias matters when the authors:

- Compute **mean climate variables per species**,
- Use those means to assign **biomes**, and
- Later compare **trait-space positions between biomes** or link trait diversity to productivity (NPP).

It doesn’t directly affect the PCA of traits (which uses trait values only), but it affects **how the PCA scores are interpreted by biome or environment** — e.g., whether “tropical” vs “temperate” species differ.

---

### **3. Why “overrepresenting well-sampled regions” matters**

Imagine two species:

- Species A lives mostly in tropical forest but is well-sampled in botanical gardens in Europe.
- Species B lives mostly in arid savannas but has thousands of records in agricultural plots in North America.

If you just average all GBIF records, both species might appear to come from *temperate* climates, even though that’s not their native biome.

This distorts the relationships between **trait syndromes** and **biome or climate**, making it look like biomes overlap more (or less) than they really do.

---

### **4. How “averaging per biome” helps**

To reduce that sampling bias, instead of pooling all occurrence records directly, you can:

1. Compute **the mean climate for each species**, weighting each biome equally (so well-sampled regions don’t dominate).
2. Or, more simply, take **the median or mean of the species-level climate means** within each biome, giving each species equal weight.

That way, a species with 100,000 records in Europe doesn’t count 100,000× more than a species with 10 records in the tropics.

---

### **5. Alternative improvements**

- **Spatial thinning**: randomly subsample GBIF occurrences so each species has an even spatial spread (e.g., one point per 50 km).
- **Sensitivity check**: test whether biome-level trait differences remain similar after thinning or equal-weight averaging.

---

### **6. In short**

> GBIF data were used to assign species to biomes and extract their climate context.
> 
> 
> But because GBIF sampling is uneven across regions, some species’ climate means (and hence biome assignments) may be biased toward well-sampled areas.
> 
> This can distort the interpretation of trait–biome relationships (even if the trait PCA itself is unaffected).
> 
> Averaging climate per species or per biome with equal weighting — or spatially thinning the GBIF data — helps ensure that well-sampled regions don’t dominate the environmental signal.
>