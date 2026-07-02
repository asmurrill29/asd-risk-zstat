# The Calculation and Validation of a TADA Posterior-Derived Z-statistic, Z'
---

## Overview

The TADA statistical model, or Transmission And De novo Association, was developed as a way of identifying risk
genes for Autism Spectrum Disorder (ASD) through a Bayesian likelihood model from WES mutation data annotated as *de novo* and inherited
mutations. TADA+ was derived to extend this model to include mutations of increasing deleteriousness (synonymous, missense, frameshift, etc.).
In both models, Bayes Factors are calculated and a subsequent FDR-like q-value per gene is used as a final metric of relative risk. My research of ASD rare-variant-enriched
genes in downstream pathway analysis (see "common-vs-rare-variants-ASD" repo) demonstrated that although this rare-variant driven analysis of risk is valuable and effective, the 
directionless, 0-1 bounded nature of its output statistic constrains the biological relevance of the TADA/TADA+ output. 

Thus, I have proposed the usage of a directed, normally distributed Wald Z-like statistic, provisionally named "Z-prime" (Z'), derived from the posterior probability of the TADA/TADA+ as mathematically supported by the 
Bernstein von Mises theorem (BvM). Here I will use develop a Stan Hamiltonian Monte Carlo algorithm to estimate the relative risk distribution for each tested gene, from which the posterior mean and standard deviation can
be calculated. From then, Z' can be generated per gene. 

A developing validation step is to compare previous q-value derived fgsea pathway analysis results with that of Z'. 

[Full theoretical derivation of Z'](docs/zprime_derivation.pdf) 

---

## Project Structure

```
project-root/
├── README.md
├── data/
│   └── raw/
│       └── satterstrom_counts.xlsx
├── scripts/
│   └── zprime_calc.R         # Generation of posterior distribution/z-prime calculation
└── results/
```

---

## Data Sources

| Dataset | Access |
|-------------|--------|
| Satterstrom et al. 2020 ASD-tested genes with their mutation counts/rates | [Supplementary Tables](https://pmc.ncbi.nlm.nih.gov/articles/PMC7250485/#SD6) |


---

## Dependencies

```r
# Install required packages
install.packages("rstan")
```
---

## Limitations

- Currently performed only on *de novo* protein truncating variants (dn-PTVs)
- Z' assumes sufficient sample size for BvM approximation to hold

---

## References

1. He, X., Sanders, et al. (2013). Integrated Model of De Novo and Inherited Genetic Variants Yields Greater Power to Identify Risk Genes. PLoS Genetics, 9(8), e1003671. https://doi.org/10.1371/journal.pgen.1003671
2. Satterstrom, F. K. et al. (2020). Large-Scale Exome Sequencing Study Implicates Both Developmental and Functional Changes in the Neurobiology of Autism. *Cell* 180(3), 568–584.

---
