# Efficient Inference for Incremental Causal Effects of Time to Treatment

Code accompanying **["Efficient Inference for Incremental Causal Effects of Time
to Treatment"](https://arxiv.org/abs/2605.29348)** (**Zhao, Z.**, Ying, A.,  and Xu, R.; under review at *Biometrika*).

This paper builds on the incremental-causal-effect framework for continuous time
to treatment and derives the **efficient influence function (EIF)** for the
estimand $\psi(\theta) = E[Y_{T(\theta)}]$, where the treatment-initiation
hazard is scaled by a factor $\theta$. The EIF is Neyman-orthogonal, so the
resulting AIPW-type estimator admits $\sqrt{n}$-inference **even when the
nuisance functions are estimated by flexible machine learning**, via
cross-fitting. The method also yields **uniform confidence bands** over a range
of $\theta$ using a multiplier bootstrap.

## Estimator

The AIPW estimator averages the influence function $\phi$ over the sample,

$$\hat{\psi}(\theta) = \frac{1}{n} \sum_{i=1}^{n} \phi_i\big(\theta; \hat{\Lambda}, \hat{\mu}\big),$$

with two nuisances:

- $\Lambda(t|l)$ — the treatment-initiation cumulative hazard,
- $\mu(u, l) = E[Y | U = u, L = l]$ — the outcome regression.

Because the score is Neyman-orthogonal, ML nuisance estimates combined with
cross-fitting preserve root-$n$ inference.

## Repository structure

```
incremental-causal-eif/
├── README.md
├── R/
    ├── simulation_eif.R        # Section 5: AIPW with (semi)parametric nuisances, NO cross-fitting
    ├── simulation_ml.R         # Section 5: AIPW with ML nuisances + K-fold cross-fitting
    └── data_application_ml.R   # Section 6: HPV / CIN2+ application (cross-fitted ML estimator)
```

## The three scripts

### `R/simulation_eif.R` — semiparametric nuisances, no cross-fitting
Evaluates the EIF/AIPW estimator $\hat{\psi}(\theta)$ with **Cox** for the
hazard $\Lambda$ and a **linear model** for $\mu$. Pointwise inference. Run per
sample size / seed:

```bash
Rscript simulation_eif.R <N> <seed>     # writes eif_<N>/seed_<seed>.rds
```

### `R/simulation_ml.R` — ML nuisances with cross-fitting
Same EIF, but the nuisances are estimated flexibly — **spline hazard regression
(`hare`)** for $\Lambda$ and a **random forest (`ranger`)** for $\mu$ — and
combined with **K-fold cross-fitting** ($K = 10$). This is the setting that lets
ML nuisances keep $\sqrt{n}$-inference. Pointwise inference. Run:

```bash
Rscript simulation_ml.R <N> <seed>      # writes ml_<N>/seed_<seed>.rds
```

Reproducing the paper's tables means running each simulation script across many
seeds (R = 1000 replications) and sample sizes $n \in \{200, 1000, 5000\}$ and
averaging the per-seed `.rds` files.

### `R/data_application_ml.R` — HPV / CIN2+ application (Section 6)
Applies the **cross-fitted ML estimator** to the cervical-cancer-screening
example, estimating $\psi(\theta)$ for the PreTectProofer group over a grid of
$\theta$, with
**pointwise 95% CIs and a uniform 95% confidence band** (multiplier bootstrap),
and produces the comparison figure.

## Data

The original Norwegian screening data are not publicly available. Following
Section 6 of the paper, the application uses the **simulated** data set from the
companion repository of Røysland et al. (2025):

> https://github.com/palryalen/paper-code

Download their simulated data, save it as `sim_data.RData` in the working
directory, and then run `data_application_ml.R`. The simulated data has 1428
individuals (1147 Amplicor/HC2, 281 PreTectProofer). No data file is included in
this repository.

## Requirements

R with:

```r
install.packages(c("survival", "flexsurv", "polspline", "ranger",
                   "dplyr", "ggplot2", "latex2exp"))
```

`simulation_eif.R` uses `survival`, `flexsurv`; `simulation_ml.R` uses
`polspline`, `ranger`; `data_application_ml.R` uses all of the above plus
`dplyr`, `ggplot2`, `latex2exp`.

## Related

The IPW method this paper improves on (Ying, Zhao & Xu, ICLR 2025) is in a
separate repository.

