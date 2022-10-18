---
name: MixedModels
topic: Mixed, Multilevel, and Hierarchical Models in R
maintainer: Ben Bolker, Julia Piaskowski, Emi Tanaka, Phillip Alday, Wolfgang Viechtbauer
e-mail: bolker@mcmaster.ca
version: 2022-10-18
source: https://github.com/cran-task-views/MixedModels/
---

**Contributors**: Ben Bolker, Michael Agronah, Julia Piaskowski, Emi Tanaka, Phillip Alday, Wolfgang Viechtbauer.


*Mixed* (or *mixed-effect*) *models* are a broad class of statistical models used to analyze data where observations can be assigned *a priori* to discrete groups, and where the parameters describing the differences between groups are treated as random (or *latent*) variables. They are one category of *multilevel*, or *hierarchical* models; *longitudinal* data are often analyzed in this framework. In econometrics, longitudinal or cross-sectional time series data are often referred to as *panel data* and are sometimes fitted with mixed models.  Mixed models can be fitted in either frequentist or Bayesian frameworks.

This task view only includes models that incorporate *continuous* (usually although not always Gaussian) latent variables. This excludes packages that handle hidden Markov models, latent Markov models, and finite (discrete) mixture models (some of these are covered by the `r view("Cluster")` task view). Dynamic linear models and other state-space models that do not incorporate a discrete grouping variable are also excluded (some of these are covered by the `r view("TimeSeries")` task view). Bioinformatic applications of [mixed models hosted on Bioconductor](https://bioconductor.org/help/search/index.html?q="mixed+models"/) are excluded as well.


### Basic model fitting

#### Linear mixed models

Linear mixed models (LMM) make the following assumptions:

- The expected values of the responses are linear combinations of the fixed predictor variables and the random effects.
- The conditional distribution of the responses is Gaussian (equivalently, the errors are Gaussian).
- The random effects are normally distributed.

*Frequentist:*

The most commonly used packages and/or functions for frequentist LMMs are:

- `r pkg("nlme", priority = "core")`: `nlme::lme()` provides REML or ML estimation. Allows multiple nested random effects, and provides structures for modeling heteroscedastic and/or correlated errors. Wald estimates of parameter uncertainty.
- `r pkg("lme4", priority = "core")`: `lmer4::lmer()`) provides REML or ML estimation. Allows multiple nested or crossed random effects, can compute profile confidence intervals and conduct parametric bootstrapping.
- `r pkg("mbest")`: fits large nested LMMs using a fast moment-based approach.

*Bayesian:*

Most Bayesian R packages use Markov chain Monte Carlo (MCMC) estimation: `r pkg("MCMCglmm", priority = "core")`, `r pkg("rstanarm")`, and `r pkg("brms", priority = "core")`; the latter two packages use the [Stan](mc-stan.org) infrastructure. `r pkg("blme")`, built on `r pkg("lme4", priority = "core")`, uses maximum a posteriori (MAP) estimation. `r pkg("bamlss")` provides a flexible set of modular functions for Bayesian regression modelling using the JAGS infrastructure.

#### Generalized linear mixed models

Generalized linear mixed models (GLMMs) can be described as hierarchical extensions of generalized linear models (GLMs), or as extensions of LMMs to different response distributions, typically in the exponential family. The random-effect distributions are typically assumed to be Gaussian on the scale of the linear predictor.

*Frequentist:*

- `r pkg("MASS")`: `MASS::glmmPQL()` fits via penalized quasi-likelihood.
- `r pkg("lme4", priority = "core")`: `lme4::glmer()` uses Laplace approximation and adaptive Gauss-Hermite quadrature; fits negative binomial as well as exponential-family models.
- `r pkg("glmmTMB", priority = "core")` uses Laplace approximation; allows some correlation structures; fits some non-exponential families (Beta, COM-Poisson, etc.) and zero-inflated/hurdle models.
- `r pkg("GLMMadaptive")` uses adaptive Gauss-Hermite quadrature; fits exponential family, negative binomial, beta, zero-inflated/hurdle/censored Gaussian models, user-specified log-densities.
- `r pkg("hglm")` fits hierarchical GLMs using $h$-likelihood (*sensu* Lee and Nelder).
- `r pkg("glmm")` fits GLMMs using Monte Carlo likelihood approximation.
- `r pkg("glmmEP")` fits probit mixed models for binary data by expectation propagation.
- `r pkg("mbest")`: fits large nested GLMMs using a fast moment-based approach.

*Bayesian:*

Most Bayesian mixed model packages use some form of Markov chain Monte Carlo (or other Monte Carlo methods).

- `r pkg("MCMCglmm", priority = "core")`: Gibbs sampling. Exponential family, multinomial, ordinal, zero-inflated/altered/hurdle, censored, multimembership, multi-response models. Pedigree (animal/kinship/phylogenetic) models.
- `r pkg("rstanarm")` Hamiltonian Monte Carlo (based on [Stan](http://mc-stan.org)); designed for `lme4` compatibility.
- `r pkg("brms", priority = "core")`: Hamilton Monte Carlo. Linear, robust linear, count data, survival, response times, ordinal, zero-inflated/hurdle/censored data.
- `r pkg("bamlss")`: optimization and derivative-based Metropolis-Hastings/slice sampling. Wide range of distributions and link functions.

Two packages (in addition to `r pkg("bamlss")`) find maximum *a posteriori* fits to Bayesian (G)LMMs by optimization:

- `r pkg("blme")` wraps `r pkg("lme4", priority = "core")` to add prior distributions.
- `r github("inbo/INLA")` provides a wide range of latent models (especially for spatial estimation), priors, and distributions.

#### Nonlinear mixed models

Nonlinear mixed models incorporate arbitrary nonlinear responses that cannot be accommodated in the framework of GLMMs. Only a few packages can accommodate *generalized* nonlinear mixed models (i.e., parametric nonlinear mixed models with non-Gaussian responses). However, many packages allow smooth nonparametric components (see ["Additive models"](#additive-models) below). Otherwise, users may need to implement GNLMMs themselves in a more general [hierarchical modeling framework](#hierarchical-modeling-frameworks).

*Frequentist:*

- `nlme::nlme()` from `r pkg("nlme")` and `lmer4::nlmer()` from `r pkg("lme4", priority = "core")` fit nonlinear mixed effects models by maximum likelihood.
- `gnlmm()` and `gnlmm3()` from `r pkg("repeated")` fit GNLMMs by Gauss-Hermite integration.
- `r pkg("saemix")` uses a stochastic approximation of the EM algorithm to fit a wide range of GNLMMs.

*Bayesian:*

- `r pkg("brms")` supports GNLMMs.

#### Generalized estimating equations

General estimating equations (GEEs) are an alternative approach to fitting clustered, longitudinal, or otherwise correlated data. These models produce estimates of the *marginal* effects (averaged across the group-level variation) rather than *conditional* effects (conditioned on group-level information).

- `r pkg("geepack", priority = "core")`, `r pkg("gee")`, and `r pkg("geeM")` are standard GEE solvers, providing GEE estimation of the parameters in mean structures with possible correlation between the outcomes.
- `r pkg("wgeesel")` implements a weighted extension of generalized linear models for longitudinal clustered data by incorporating the correlation within-cluster when data is missing at random.
- `r pkg("geesmv")`: GEE estimator using the original sandwich variance estimator proposed by Liang and Zeger (1986), and eight types of variance estimators for improving the finite small-sample performance.
- `r pkg("multgee")` is a GEE solver for correlated nominal or ordinal multinomial responses.

### Specialized models

- [**Additive models**]{#additive-models} (models incorporating smooth functional components such as regression splines or Gaussian processes): `r pkg("gamm4")`, `r pkg("mgcv")`, `r pkg("brms", priority = "core")`, `r pkg("lmeSplines")`, `r pkg("bamlss")`, `r pkg("gamlss")`, `r github("Biometris/LMMsolver")`, `r pkg("R2BayesX")`, `r pkg("GLMMRR")`.
- **Big data/distributed computation**: `r pkg("lmmpar")`, `r pkg("mbest")`. See also [MixedModels.jl](https://juliastats.org/MixedModels.jl/dev/) (Julia), [diamond](https://github.com/stitchfix/diamond) (Python).
- **Censored data** (response data known only up to lower/upper bounds): `r pkg("brms", priority = "core")` (general), `r pkg("ARpLMEC")` (censored Gaussian, autoregressive errors). Censored Gaussian (Tobit) responses: `r pkg("GLMMadaptive")`, `r pkg("MCMCglmm", priority = "core")`, `r pkg("gamlss")`.
- **Differential equations** (fitting DEs with group-structured parameters): `r pkg("mixedsde")`, see also `r view("DifferentialEquations")`.
- **Factor analytic, latent variable, and structural equation models**:  `r pkg("lavaan", priority = "core")`, `r pkg("nlmm")`,`r pkg("sem")`, `r pkg("piecewiseSEM")`, `r pkg("semtree")`, and  `r pkg("blavaan")`; also see the `r view("Psychometrics")` task view.
- **Kinship-augmented models** (responses where individuals have a known family relationship): `r pkg("pedigreemm")`, `r pkg("coxme")`, `r pkg("kinship2")`, `r github("Biometris/LMMsolver")`, `r pkg("MCMCglmm", priority = "core")`.
- **Location-scale models**: `r pkg("nlme", priority = "core")`, `r pkg("glmmTMB", priority = "core")`, `r pkg("brms", priority = "core")`, `r pkg("mgcv")` [with `family` chosen from one of the `*ls`/`*lss` options]  all allow modeling of the dispersion/scale component.
- **Missing values**: `r pkg("mice")`, `r pkg("mlmmm")` (EM imputation), `r pkg("CRTgeeDR")`, `r pkg("JointAI")`, `r pkg("mdmb")`, `r pkg("pan")`; also see the `r view("MissingData")` task view.
- **Multimembership models**: (Bayesian) `r pkg("MCMCglmm", priority = "core")`, `r pkg("brms", priority = "core")`, `r github("benrosche/rmm")`; (frequentist) `r github("jvparidon/lmerMultiMember")`.
- **Multinomial responses**: `r pkg("bamlss")`, `r pkg("R2BayesX")`, `r pkg("MCMCglmm", priority = "core")`, `r pkg("mgcv")`.
- **Multi-trait analysis**: (multiple dependent variables) `r pkg("BMTME")`, `r pkg("MCMCglmm", priority = "core")`.
- **Non-Gaussian random effects**: `r pkg("brms", priority = "core")`, `r pkg("repeated")`, `r pkg("spaMM")`.
- **Ordinal-valued responses** (responses measured on an ordinal scale): `r pkg("ordinal")`, `r pkg("cplm")`.
- **Over-dispersed models**: `r pkg("aod")`, `r pkg("aods3")`.
- **Panel data**: in econometrics, *panel data* typically refers to subjects (individuals or firms) that are sampled repeatedly over time. The theoretical and computational approaches used by econometricians overlap with mixed models (e.g., see [here](https://cran.r-project.org/web/packages/plm/vignettes/A_plmPackage.html#nlme)). The `r pkg("plm")` package can fit mixed-effects panel models; see also the `r view("Econometrics")` task view.
- **Quantile regression**: `r pkg("lqmm")`, `r pkg("qrNLMM")`.
- **Phylogenetic models**: `r pkg("pez")`, `r pkg("phyr")`, `r pkg("MCMCglmm", priority = "core")`, `r pkg("brms", priority = "core")`.
- **Regularized/penalized models** (regularization or variable selection by ridge, lasso, or elastic net penalties): `r pkg("splmm")` fits LMMs for high-dimensional data by imposing penalty on both the fixed effects and random effects for variable selection. `r pkg("glmmLasso")` fits GLMMs with L1-penalized (LASSO) fixed effects. `r pkg("bamlss")` implements LASSO-like penalization for generalized additive models.
- **Robust/heavy-tailed estimation** (downweighting the importance of extreme observations): `r pkg("robustlmm")`, `r pkg("robustBLME")` (Bayesian robust LME), `r pkg("CRTgeeDR")` for the doubly robust inverse probability weighted augmented GEE estimator. Some packages (`r pkg("brms", priority = "core")`, `r pkg("bamlss")`, `r pkg("mgcv")` with `family = "scat"`) allow heavy-tailed response distributions such as Student-$t$.
- **Skewed data**: `r pkg("skewlmm")` fits a scale mixture of skew-normal linear mixed models using expectation-maximization (EM).
- **Spatial models**: `r github("inbo/INLA")`, `r pkg("nlme", priority = "core")` (with `corStruct` functions), `r pkg("CARBayesST")`, `r pkg("sphet")`, `r pkg("spind")`, `r pkg("spaMM")`, `r pkg("glmmfields")`, `r pkg("glmmTMB")`, `r pkg("inlabru")` (spatial point processes via log-Gaussian Cox processes), `r pkg("brms", priority = "core")`, `r github("Biometris/LMMsolver")`, `r pkg("bamlss")`; see also the `r view("Spatial")` and `r view("SpatioTemporal")` CRAN task views.
- **Survival analysis**: `r pkg("coxme")`.
- **Tree-based models**: `r pkg("glmertree")`, `r pkg("semtree")`.
- **Zero-inflated models**: (frequentist) `r pkg("glmmTMB")`, `r pkg("cplm")`; (Bayesian): `r pkg("MCMCglmm", priority = "core")`, `r pkg("brms", priority = "core")`, `r pkg("bamlss")`, `r pkg("mgcv")` (zi Poisson only).
- **Bioinformatics/quantitative genetics**: `r pkg("MCMC.qpcr")`, `r pkg("QGglmm")`, `r pkg("CpGassoc")` (methylation studies).
- **By area**:
    * `r pkg("mvglmmRank")`, multivariate generalized linear mixed models for ranking sports teams.


### Hierarchical modeling frameworks

These packages do not directly provide functions to fit mixed models, but instead implement interfaces to general-purpose sampling and optimization toolboxes that can be used to fit mixed models. While models require extra effort to set up, and often require programming in a domain-specific language other than R, these frameworks are more flexible than most of the other packages listed here.

* Interfaces to [JAGS](https://mcmc-jags.sourceforge.io/)/[OpenBUGS](https://www.mrc-bsu.cam.ac.uk/software/bugs/openbugs/): `r pkg("R2jags")`, `r pkg("rjags")`, `r pkg("R2OpenBUGS")` (BUGS language).
* Interfaces to [Stan](http://mc-stan.org/) (C++ extensions): `r pkg("rstan")`, `r github("stan-dev/cmdstanr")`, `r github("rmcelreath/rethinking")`.
* Other frameworks: `r pkg("TMB")` (automatic differentiation and Laplace approximation) (C++ extensions), `r pkg("tmbstan")`, `r pkg("nimble")`, `r pkg("greta")` (R interface to TensorFlow).


### Model diagnostics and summary statistics

#### Model diagnostics

- **general**: `r pkg("HLMdiag")` (diagnostic tools for hierarchical (multilevel) linear models), `r pkg("rockchalk")`, `r pkg("performance")`, `r pkg("multilevelTools")`, `r pkg("merTools")` (for models fitted using `lme4`).
- **influential data points**: `r pkg("influence.ME")`, `r pkg("influence.SEM")`.
- **residuals**: `r pkg("DHARMa")`.

#### Summary statistics

- **Correlations**:  `r pkg("iccbeta")` (intraclass correlation), `r pkg("r2glmm")` (R^2 and partial R^2).
- **Information criteria**: `r pkg("cAIC4")` (conditional AIC) , `r pkg("blmeco")` (WAIC).
- **Robust variance-covariance estimates**: `r pkg("clubSandwich")`, `r pkg("merDeriv")`.

#### Derivatives

The first and second derivatives of log-likelihood with respect to parameters can be useful for various model evaluation tasks (e.g., computing sensitivities, robust variance-covariance matrices, or delta-method variances).

- `r pkg("lmeInfo")`, `r pkg("merDeriv")`.


### Data sets

Many packages include small example data sets (e.g., `r pkg("lme4", priority = "core")`, `r pkg("nlme", priority = "core")`). These packages provide previously described data sets often used in evaluating mixed models.

- `r pkg("mlmRev")`: examples from the Multilevel Software Comparative Reviews.
- `r pkg("SASmixed")`: data sets from *[SAS System for Mixed Models](https://v8doc.sas.com/sashtml/hrddoc/indfiles/55235.htm)*.
- `r pkg("StroupGLMM")`: R scripts and data sets for *[Generalized Linear Mixed Models](https://www.taylorfrancis.com/books/mono/10.1201/b13151/generalized-linear-mixed-models-walter-stroup)*.
- `r pkg("blmeco")`: Data and functions accompanying *[Bayesian Data Analysis in Ecology using R, BUGS and Stan](https://www.elsevier.com/books/bayesian-data-analysis-in-ecology-using-linear-models-with-r-bugs-and-stan/korner-nievergelt/978-0-12-801370-0)*.
- `r pkg("nlmeU")`: Data sets, functions and scripts described in *[Linear Mixed-Effects Models: A Step-by-Step Approach](https://link.springer.com/book/10.1007/978-1-4614-3900-4)*.
- `r pkg("VetResearchLMM")`: R scripts and data sets for *[Linear Mixed Models. An Introduction with applications in Veterinary Research](https://hdl.handle.net/10568/5379)*.
- `r pkg("languageR")`: R scripts and data sets for *[Analyzing Linguistic Data: A practical introduction to statistics using R](https://doi.org/10.1017/CBO9780511801686)*.


### Model presentation and prediction

Functions and frameworks for convenient and tabular and graphical output of mixed model results:

- **Tables**: `r pkg("huxtable")`, `r pkg("broom.mixed", priority = "core")`, `r pkg("rockchalk")`, `r pkg("parameters")`, `r pkg("modelsummary")`.
- **Figures**: `r pkg("dotwhisker")`, `r pkg("sjPlot")`, `r pkg("rockchalk")`.


### Convenience wrappers

These functions provide convenient frameworks to fit and interpret mixed models.

- **Model fitting**: `r pkg("multilevelmod", priority = "core")`,  `r pkg("ez")`, `r pkg("mixlm")`, `r pkg("afex")`, `r pkg("dalmatian")` (wrapper to JAGS and `r pkg("nimble")`).
- **Model summary**: `r pkg("broom.mixed", priority = "core")`, `r pkg("insight")`.
- **Variable selection & model averaging**: `r pkg("LMERConvenienceFunctions")`, `r pkg("MuMIn")`, `r pkg("glmulti")` (see, e.g., [maintainer's blog](https://vcalcagnoresearch.wordpress.com/package-glmulti/) or [here](https://gist.github.com/bbolker/4ae3496c0ddf99ea2009a22b94aecbe5) for use with mixed models).


### Inference and model selection

#### Hypothesis testing

- **Fixed effects**: `r pkg("car")`, `r pkg("lmerTest")`, `r pkg("RVAideMemoire")`, `r pkg("emmeans")`, `r pkg("afex")`, `r pkg("pbkrtest")`, `r pkg("CLME")`.
- **Random effects**: `r pkg("varTestnlme")`, `r pkg("RLRsim")`, `r pkg("mvctm")`.

#### Prediction and estimation

- `r pkg("emmeans")`, `r pkg("effects")`, `r pkg("margins")`, `r pkg("MarginalMediation")`, `r pkg("marginaleffects")`.

#### Bootstrapping

- `r pkg("pbkrtest")`, `r pkg("lme4", priority = "core")` (`lme4::bootMer()` function), `r pkg("lmeresampler")`.

#### Power analysis

- `r pkg("longpower")`, `r pkg("clusterPower")`, `r pkg("pass.lme")`, `r pkg("simr")`.

#### Model selection

- `r pkg("cAIC4")` (`cAIC4::stepcAIC`), `r pkg("buildmer")`, `r pkg("MuMIn")`, `r github("timnewbold/StatisticalModels")` (`GLMERSelect`).


### Other

**Commercial software interfaces:**

- [Mplus](https://www.statmodel.com/): `r pkg("MplusAutomation")`.
- [ASReml R](https://vsni.co.uk/software/asreml-r): `r pkg("asremlPlus")`.
- [Phoenix NLME software](http://www.certara.com/software/pkpd-modeling-and-simulation/phoenix-nlme/): `r pkg("Phxnlme")`.


### Links
- Help: [R-SIG-mixed-models mailing list](https://stat.ethz.ch/mailman/listinfo/r-sig-mixed-models) for discussion of mixed-model-related questions, course announcements, etc..
- Help: [[r] + [mixed-models] tags on Stack Overflow](http://stackoverflow.com/questions/tagged/r+mixed-models).
- Help: [Cross Validated](http://stats.stackexchange.com/).
- Other software: [Mixed models Bioconductor](https://bioconductor.org/help/search/index.html?q="mixed+models"/)
- Other software: [ASReml-R](https://vsni.co.uk/software/asreml-r) (`r pkg("asremlPlus")`).
- Other software: [assist](https://yuedong.faculty.pstat.ucsb.edu/software.html).
- Other software: [INLA](http://www.r-inla.org/home).
- Other software: [Zelig Project](https://zeligproject.org/).
- Other software: [MixWild/MixRegLS](https://voices.uchicago.edu/hedeker/mixwild_mixregls/) for scale-location modeling.
- Other software: [MixedModels.jl](https://github.com/JuliaStats/MixedModels.jl) for mixed models in Julia.
- Book: *[Mixed-Effects Models in S and S-PLUS](https://link.springer.com/book/10.1007/b98882)*.
- Book: *[SAS System for Mixed Models](https://v8doc.sas.com/sashtml/hrddoc/indfiles/55235.htm)*.
- Book: *[Generalized Linear Mixed Models](https://www.taylorfrancis.com/books/mono/10.1201/b13151/generalized-linear-mixed-models-walter-stroup)*.
- Book: *[Bayesian Data Analysis in Ecology using R, BUGS and Stan](https://www.elsevier.com/books/bayesian-data-analysis-in-ecology-using-linear-models-with-r-bugs-and-stan/korner-nievergelt/978-0-12-801370-0)*.
- Book: *[Linear Mixed-Effects Models: A Step-by-Step Approach](https://link.springer.com/book/10.1007/978-1-4614-3900-4)*.
- Book: *[Mixed Effects Models and Extensions in Ecology with R](https://link.springer.com/book/10.1007/978-0-387-87458-6)*.
- Online Book: *[Embrace Uncertainty: Mixed-effects models with Julia](https://juliamixedmodels.github.io/EmbraceUncertainty/)*.
