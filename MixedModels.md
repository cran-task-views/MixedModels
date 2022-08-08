---
name: MixedModels
topic: Mixed, multilevel, and hierarchical models in R
maintainer: Ben Bolker, Julia Piaskowski
e-mail: bolker@mcmaster.ca
version: 2022-07-29
source: https://github.com/bbolker/mixedmodels-misc/blob/master/taskview/MixedModels.md
---

**Authors**: Ben Bolker, Michael Agronah, ??  

*Mixed models* are a broad class of statistical models used to analyze data where observations can be assigned to discrete groups, and where the parameters describing the differences are treated as *random variables*. They are also described as *multilevel*, or *hierarchical*,  models; *longitudinal* data are often analyzed in this framework.  Mixed models can be fitted in either frequentist or Bayesian frameworks.

**Scope**: only including models that incorporate *continuous* (usually although not always Gaussian) latent variables; this excludes packages that handle hidden Markov Models, finite (discrete) mixture models, latent Markov models, etc.

## Basic model fitting

This section is to describe linear mixed models (LMM) under the following assumpitions:
    1. The responses are linear combinations of the predictor variables;
    2. the error residuals are normally-distributed residuals;
    3. the random effects are normally distributed.  

### Linear mixed models

#### Frequentist

The most commonly used packages and/or functions for frequentist LMMs are:

- `r pkg("nlme", priority = "core")`: `nlme::lme` conducts REML or ML estimation, can include multiple nested random effects, and construction residual correlation structures for heteroscedasticity.
- `r pkg("lme4", priority = "core")`: `lmer4::lmer`) conducts REML or ML estimation, can included multiple nested and random effects, can profile confidence intervals, and conduct parametric bootstrapping. 

#### Bayesian

Most Bayesian R packages use Markov chain Monte Carlo (MCMC) estimation: `r pkg("MCMCglmm", priority = "core")`, `r pkg("rstanarm")` and `r pkg("brms", priority = "core")`; the latter two packages uses the [Stan](mc-stan.org) infrastructure. `r pkg("blme")`, built on `r pkg("lme4")`, uses maximum a posteriori (MAP) estimation.

### Generalized linear mixed models

Generalized linear mixed models (GLMMs) can be described as hierarchical extensions of generalized linear models (GLMs), or a extensions of LMMs to different response distributions, typically in the exponential family. The random-effect distributions are typically assumed to be Gaussian on the scale of the linear predictor.

#### Frequentist

- `r pkg("MASS")`: `MASS::glmmPQL` fits via penalized quasi-likelihood.
- `r pkg("lme4")`: `lme4::glmer` does a Laplace approximation and adaptive Gauss-Hermite quadrature, and `r pkg("glmmTMB")` also does a Laplace approximation).
- `r pkg("GLMMadaptive")` and `r pkg("hglm")` handle hierarchical GLMs.
- `r pkg("lmeNB")` implements a negative binomial distribution
  
#### Bayesian
 
- `r pkg("MCMCglmm")` fits GLMMs using MCMC techniques. 
- `r pkg("rstanarm")` fits GLMMs using Markov Chain Monte Carlo, variational approximations to the posterior distribution, or optimization. 
- `r pkg("brms")` supports a wide range of distributions and link functions for fitting GLMMs. 
- `r pkg("glmm")` fits GLMMs using Monte Carlo Likelihood Approximation.
- *binary data* Two packages can handle binary data. `r pkg("glmmEP")`, which handles probit models and `r pkg("GLMMRR")` which can use one of four different cumulative distribution functions. 

### Nonlinear mixed models

Nonlinear mixed models incorporate arbitrary nonlinear responses that cannot be accommodated in the framework of GLMMs. Only a few packages can accommodate **generalized** nonlinear mixed models
(i.e., nonlinear mixed models with non-Gaussian responses).

#### Frequentist

- The functions `nlme::nlme()` from `r pkg("nlme")`, `lmer4::nlmer()` from `r pkg("lme4")` and `GNLMM()` from `r pkg("repeated")` can conduct basic model fitting.
- `r pkg("saemix")` provides a stochastic approximation of the EM algorithm.

#### Bayesian

- `r pkg("brms")` supports non-linear mixed models. 

### Generalized estimating equations 

General estimating equations (GEEs) are an alternative approach to fitting clustered, longitudinal, or otherwise correlated data. These models produce estimates of the *marginal* effects (averaged across the group-level variation) rather than *conditional* effects (conditioned on group-level information).

- `r pkg("geepack", priority = "core")`, `r pkg("gee")` and `r pkg("geeM")` are standard GEE solvers, providing GEE estimation of the parameters in mean structures with possible correlation between the outcomes. 
- `r pkg("wgeesel")` implements a weighted extension of generalized linear models for longitudinal clustered data by incorporating the correlation within-cluster when data is missing at random. The parameters in mean, scale correlation structures are estimated based on quasi-likelihood. 
- `r pkg("geesmv")`: GEE estimator using the original sandwich variance estimator proposed by Liang and Zeger ([1986](http://ibg.colorado.edu/cdrom2011/medland/fri2011/HWSE.pdf)), and eight types of variance estimators for improving the finite small-sample performance. 
- `r pkg("multgee")` is a GEE solver for correlated nominal or ordinal multinomial responses using a local odds ratios parameterization.

## Specialized models

- **Censored data**: `r pkg("brms")` (general), `r pkg("lmec")` (censored Gaussian), `r pkg("ARpLMEC")` (censored Gaussian, autoregressive errors) and `r pkg("tlmec")` (censored Gaussian and Student-t distributions)

- **Differential equations**: `r pkg("mixedsde")`, `r pkg("nlmeODE")` `r pkg("PSM")`, see also `r view("DifferentialEquations")`

- **Factor analytic, latent variable, and structural equation models**:  `r pkg("lavaan", priority = "core")`, `r pkg("nlmm")`,`r pkg("sem")`, `r pkg("piecewiseSEM")`, `r pkg("semtree")`, `r pkg("semPLS")` and  `r pkg("blavaan")`. (See also the `r view("Psychometrics")` task view)

- **kinship-augemented models**: `r pkg("pedigreemm")`, `r pkg("coxme")`, `r pkg("kinship2")`

- **Missing values**: `r pkg("mlmmm")` (EM imputation), `r pkg("CRTgeeDR")`, also see the `r view("MissingData")` task view for strategies for imputing missing data

- **Multinomial responses**: FIXME

-- **multi-trait analysis**: (multiple dependent variables) `r pkg("BMTME")`

- **Ordinal-valued responses**: `r pkg("ordinal")`, `r pkg("cplm")`

- **Quantile regression**: `r pkg("lqmm")`, `r pkg("qrLMM")`,`r pkg("qrNLMM")`

- **Phylogenetic linear mixed models**: `r pkg("pez")`

- **Regularized/Penalized models** (regularization or variable selection by ridge, lasso, or elastic net penalties): `r pkg("splmm")` fits LMMs for high-dimensional data by imposing penalty on both the fixed effects and random effects for variable selection.

- **Robust estimation** for downweighting the importance of extreme observations: `r pkg("robustlmm")`, `r pkg("robustBLME")` (Bayesian robust LME), `r pkg("CRTgeeDR")` for the doubly robust inverse probability weighted augmented GEE estimator. 

- **Survival analysis**: `r pkg("coxme")`

- **Spatial models**: `r github("inbo/INLA")`, `r pkg("nlme")` (with `corStruct` functions), `r pkg("CARBayesST")`, `r pkg("sphet")`, `r pkg("spind")`, `r pkg("spaMM")`, `r pkg("glmmfields")`, `r pkg("glmmTMB")`, `r pkg("inlabru")` (spatial point processes via log-Gaussian Cox processes), `r pkg("brms")`; also see the `r view("Spatial")` and `r view("SpatioTemporal")` CRAN task views

- **Tree-based models**: `r pkg("glmertree")`, `r pkg("semtree")`

- **Zero-inflated models**: (frequentist) `r pkg("glmmTMB")`, `r pkg("cplm")`; (Bayesian): `r pkg("MCMCglmm")`, `r pkg("brms")`

## Model diagnostics and summary statistics

### Model diagnostics

`r pkg("HLMdiag")`, `r pkg("rockchalk")`, `r pkg("influence.ME")`, `r pkg("aods3")` (overdispersion), `r pkg("DHARMa")`, `r pkg("performance")`

### Summary statistics

`r pkg("iccbeta")` (intraclass correlation), `r pkg("r2glmm")` (R^2 and partial R^2),
`r pkg("HiLMM")` (heritability), `r pkg("cAIC4")` (conditional AIC) , `r pkg("blmeco")` (WAIC)

### Derivatives

`r pkg("lmeInfo")`, `r pkg("merDeriv")`, `r pkg("lmmpar")`

**robust variance-covariance estimates**: `r pkg("clubSandwich")`, `r pkg("merDeriv")`

## Datasets

`r pkg("mlmRev")`, `r pkg("lme4")`, `r pkg("nlme")`, `r pkg("SASmixed")`, `r pkg("StroupGLMM")`, `r pkg("blmeco")`, `r pkg("nlmeU")`, `r pkg("VetResearchLMM")`

## Model presentation and prediction

Functions and frameworks for convenient and tabular and graphical output of mixed model results: `r pkg("effects")` `r pkg("emmeans")`, `r pkg("dotwhisker")`, `r pkg("huxtable")`, `r pkg("sjPlot")`, `r pkg("rockchalk")`

`r pkg("broom.mixed")`, `r pkg("insight")`

## Convenience wrappers

These functions don't necessarily add new functionality, but
provide convenient frameworks for less experienced users to fit and interpret mixed models.

 `r pkg("ez")`, `r pkg("mixlm")` `r pkg("afex")`, `r pkg("RVAideMemoire")`,
`r pkg("ZeligMultilevel")` `r pkg("cubature")`.

### Model and variable selection

`r pkg("LMERConvenienceFunctions")`, `r pkg("MuMIn")`

## Inference

`r pkg("pbkrtest")`, `r pkg("afex")` `r pkg("varTestnlme")`, `r pkg("lmeVarComp")`, `r pkg("RLRsim")`, `r pkg("car")` (`Anova()`), `r pkg("CLME")`, `r pkg("lmerTest")`

## Bootstrapping

`r pkg("pbkrtest")`, `r pkg("lme4")` (`bootMer` function), `r pkg("lmeresampler")`, `r pkg("glmmboot")`

#### Additive models

`r pkg("gamm4")`, `r pkg("mgcv")`, `r pkg("brms")`, `r pkg("lmeSplines")`

(note `assist` package is currently archived on CRAN)

### Bioinformatic applications

(FIXME: refer to/check bioconductor?)

`r pkg("MCMC.qpcr")`,`r pkg("CpGassoc")`, `r pkg("QGglmm")`, `r pkg("Phxnlme")`, `r pkg("mlmm.gwas")`

#### Power analysis

`r pkg("longpower")`, `r pkg("clusterPower")`, `r pkg("powerlmm")`, `r pkg("pass.lme")`.

### Off-CRAN, commercial tools, etc

- **interfaces to commercial software**: `r pkg("MplusAutomation")` ([Mplus](https://www.statmodel.com/)), `r pkg("asremlPlus")` ([ASReml R](https://vsni.co.uk/software/asreml-r))
- **off-CRAN, open source**: [INLA](http://www.r-inla.org/home)
- **tools for general-purpose latent variable/mixed modeling**: [JAGS](https://mcmc-jags.sourceforge.io/) (via `r pkg("jags")`/`r pkg("r2jags")`), [Stan](https://mc-stan.org) (via `r pkg("rstan")`), `r pkg("nimble")`, `r pkg("TMB")`, `r pkg ("greta")`

#### Other

The following are other packages applied in mixed models. ; `r pkg("lmeNBBayes")` `r pkg("MarginalMediation")`
`r pkg("skewlmm")`  fits scale mixture of skew-normal linear mixed models using  expectation-maximization (EM)
`r pkg("mvglmmRank")` implements multivariate Generalized Linear Mixed Models for ranking sport teams

## Links

- [R-SIG-mixed-models mailing list](https://stat.ethz.ch/mailman/listinfo/r-sig-mixed-models) for discussion of mixed-model-related questions, course announcements, etc..
- [r+mixed-models tags on Stack Overflow](http://stackoverflow.com/questions/tagged/r+mixed-models)
- [Cross Validated](http://stats.stackexchange.com)
- [INLA](http://www.r-inla.org/home)
