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

There is a large number of mixed model packages for bioinformatic applications found on [Bioconductor](https://bioconductor.org/help/search/index.html?q=mixed+models/). 
## Basic model fitting

This section is to describe linear mixed models (LMM) under the following assumpitions:
    1. The responses are linear combinations of the predictor variables;
    2. the error residuals are normally-distributed residuals;
    3. the random effects are normally distributed.  
### Linear mixed models
#### Frequentist

The most commonly used packages and/or functions for frequentist LMMs are:

- `r pkg("nlme", priority = "core")`: `nlme::lme()` conducts REML or ML estimation, can include multiple nested random effects, and construction residual correlation structures for heteroscedasticity.
- `r pkg("lme4", priority = "core")`: `lmer4::lmer()`) conducts REML or ML estimation, can included multiple nested and random effects, can profile confidence intervals, and conduct parametric bootstrapping. 
- `r pkg("TMB")` is a flexible tool to implement complex random effect models  

#### Bayesian

Most Bayesian R packages use Markov chain Monte Carlo (MCMC) estimation: `r pkg("MCMCglmm", priority = "core")`, `r pkg("rstanarm")` and `r pkg("brms", priority = "core")`; the latter two packages uses the [Stan](mc-stan.org) infrastructure. `r pkg("blme")`, built on `r pkg("lme4", priority = "core")`, uses maximum a posteriori (MAP) estimation.

Other packages: `r pkg("jags")`, `r pkg("bayesmix")`, `r pkg("R2jags")`, `r pkg("greta")`

### Generalized linear mixed models

Generalized linear mixed models (GLMMs) can be described as hierarchical extensions of generalized linear models (GLMs), or a extensions of LMMs to different response distributions, typically in the exponential family. The random-effect distributions are typically assumed to be Gaussian on the scale of the linear predictor.

#### Frequentist

- `r pkg("MASS")`: `MASS::glmmPQL()` fits via penalized quasi-likelihood.
- `r pkg("lme4", priority = "core")`: `lme4::glmer()` does a Laplace approximation and adaptive Gauss-Hermite quadrature, and `r pkg("glmmTMB")` also does a Laplace approximation).
- `r pkg("GLMMadaptive")` and `r pkg("hglm")` handle hierarchical GLMs.
- `r pkg("lmeNB")` implements a negative binomial distribution
- `r pkg("mvglmmRank")`, multivariate generalized linear mixed models for ranking sports teams
#### Bayesian
 
- `r pkg("MCMCglmm", priority = "core")` fits GLMMs using MCMC techniques. 
- `r pkg("rstanarm")` fits GLMMs using Markov Chain Monte Carlo, variational approximations to the posterior distribution, or optimization. 
- `r pkg("brms", priority = "core")` supports a wide range of distributions and link functions for fitting GLMMs. 
- `r pkg("glmm")` fits GLMMs using Monte Carlo Likelihood Approximation.
- `r pkg("MCMC.qpcr")`, quantitative RT-PCR data are fit with generalized linear mixed models and a lognormal-Poisson error distribution using MCMC
- **binary data** Two packages can handle binary data. `r pkg("glmmEP")`, which handles probit models and `r pkg("GLMMRR")` which can use one of four different cumulative distribution functions. 
### Nonlinear mixed models

Nonlinear mixed models incorporate arbitrary nonlinear responses that cannot be accommodated in the framework of GLMMs. Only a few packages can accommodate **generalized** nonlinear mixed models
(i.e., nonlinear mixed models with non-Gaussian responses).
#### Frequentist

- The functions `nlme::nlme()` from `r pkg("nlme")`, `lmer4::nlmer()` from `r pkg("lme4", priority = "core")` and `GNLMM()` from `r pkg("repeated")` can conduct basic model fitting.
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

- **Additive models**: `r pkg("gamm4")`, `r pkg("mgcv")`, `r pkg("brms", priority = "core")`, `r pkg("lmeSplines")`

- **Censored data**: `r pkg("brms", priority = "core")` (general), `r pkg("lmec")` (censored Gaussian), `r pkg("ARpLMEC")` (censored Gaussian, autoregressive errors) and `r pkg("tlmec")` (censored Gaussian and Student-t distributions)

- **Differential equations**: `r pkg("mixedsde")`, `r pkg("nlmeODE")` `r pkg("PSM")`, see also `r view("DifferentialEquations")`

- **Factor analytic, latent variable, and structural equation models**:  `r pkg("lavaan", priority = "core")`, `r pkg("nlmm")`,`r pkg("sem")`, `r pkg("piecewiseSEM")`, `r pkg("semtree")`, `r pkg("semPLS")` and  `r pkg("blavaan")`. (See also the `r view("Psychometrics")` task view)

- **kinship-augemented models**: `r pkg("pedigreemm")`, `r pkg("coxme")`, `r pkg("kinship2")`

- **latent variable/mixed modeling**: `r pkg("jags")`, `r pkg("R2jags")`, `r pkg("rstan")`), `r pkg("nimble")`, `r pkg("TMB")`, `r pkg ("greta")`

- **Missing values**: `r pkg("mlmmm")` (EM imputation), `r pkg("CRTgeeDR")`, also see the `r view("MissingData")` task view for strategies for imputing missing data

- **Multinomial responses**: FIXME

-- **multi-trait analysis**: (multiple dependent variables) `r pkg("BMTME")`

- **Ordinal-valued responses**: `r pkg("ordinal")`, `r pkg("cplm")`

- **Over-dispersed models**: `r pkg("aod")`, `r pkg("aod3")`

- **Quantile regression**: `r pkg("lqmm")`, `r pkg("qrLMM")`,`r pkg("qrNLMM")`

- **Phylogenetic linear mixed models**: `r pkg("pez")`

- **Regularized/Penalized models** (regularization or variable selection by ridge, lasso, or elastic net penalties): `r pkg("splmm")` fits LMMs for high-dimensional data by imposing penalty on both the fixed effects and random effects for variable selection.

- **Robust estimation** for downweighting the importance of extreme observations: `r pkg("robustlmm")`, `r pkg("robustBLME")` (Bayesian robust LME), `r pkg("CRTgeeDR")` for the doubly robust inverse probability weighted augmented GEE estimator. 

- **Survival analysis**: `r pkg("coxme")`

- **Spatial models**: `r github("inbo/INLA")`, `r pkg("nlme", priority = "core")` (with `corStruct` functions), `r pkg("CARBayesST")`, `r pkg("sphet")`, `r pkg("spind")`, `r pkg("spaMM")`, `r pkg("glmmfields")`, `r pkg("glmmTMB")`, `r pkg("inlabru")` (spatial point processes via log-Gaussian Cox processes), `r pkg("brms", priority = "core")`; also see the `r view("Spatial")` and `r view("SpatioTemporal")` CRAN task views

- **skewed data**: `r pkg("skewlmm")` fits scale mixture of skew-normal linear mixed models using expectation-maximization (EM)

- **Tree-based models**: `r pkg("glmertree")`, `r pkg("semtree")`

- **Zero-inflated models**: (frequentist) `r pkg("glmmTMB")`, `r pkg("cplm")`; (Bayesian): `r pkg("MCMCglmm", priority = "core")`, `r pkg("brms", priority = "core")`
## Model diagnostics and summary statistics
### Model diagnostics

- **general**:  `r pkg("HLMdiag")` (Diagnostic Tools for Hierarchical (Multilevel) Linear Models), `r pkg("rockchalk")`, `r pkg("performance")`, `r pkg("multilevelTools")`
- **influential data points**: `r pkg("influence.ME")`, `r pkg("influence.SEM")`, 
- **residuals**: `r pkg("DHARMa")`
### Summary statistics

- **Correlations**:  `r pkg("iccbeta")` (intraclass correlation), `r pkg("r2glmm")` (R^2 and partial R^2),
- **Quantitative genetics parameters**:  `r pkg("HiLMM")`, `r pkg("QGglmm")` )
- **Information criteria**: `r pkg("cAIC4")` (conditional AIC) , `r pkg("blmeco")` (WAIC)
- **robust variance-covariance estimates**: `r pkg("clubSandwich")`, `r pkg("merDeriv")`
### Derivatives

- `r pkg("lmeInfo")`, `r pkg("merDeriv")`, `r pkg("lmmpar")`
## Datasets

Many packages include data sets to provide examples to test package functions with (e.g. `r pkg("lme4", priority = "core")`, `r pkg("nlme", priority = "core")`). The packages listed here are previously described data sets often used in evaluating mixed models. 

- `r pkg("mlmRev")`: examples from the Multilevel Software Comparative Reviews
- `r pkg("SASmixed")`: data sets from *[SAS System for Mixed Models](https://v8doc.sas.com/sashtml/hrddoc/indfiles/55235.htm)*
- `r pkg("StroupGLMM")`: R scripts and data sets for *[Generalized Linear Mixed Models](https://www.taylorfrancis.com/books/mono/10.1201/b13151/generalized-linear-mixed-models-walter-stroup)*
- `r pkg("blmeco")`: Data and functions accompanying *[Bayesian Data Analysis in Ecology using R, BUGS and Stan](https://www.elsevier.com/books/bayesian-data-analysis-in-ecology-using-linear-models-with-r-bugs-and-stan/korner-nievergelt/978-0-12-801370-0)* 
- `r pkg("nlmeU")`: Data sets, functions and scripts described in *[Linear Mixed-Effects Models: A Step-by-Step Approach](https://link.springer.com/book/10.1007/978-1-4614-3900-4)* 
- `r pkg("VetResearchLMM")`: R scripts and data sets for *[Linear Mixed Models. An Introduction with applications in Veterinary Research](https://www.ilri.org/publications/linear-mixed-model-introduction-applications-veterinary-research)*
## Model presentation and prediction

Functions and frameworks for convenient and tabular and graphical output of mixed model results: 

- **Tables**: `r pkg("huxtable")`, `r pkg("broom.mixed")`, `r pkg("rockchalk")`
- **Figures**: `r pkg("dotwhisker")`, `r pkg("sjPlot")`, `r pkg("CpGassoc")` (methylation studies), `r pkg("rockchalk")` 
## Convenience wrappers

These functions provide convenient frameworks to fit and interpret mixed models.

- **Model fitting**: `r pkg("multilevelmod", priority = "core")`,  `r pkg("ez")`, `r pkg("mixlm")`, `r pkg("afex")`, `r pkg("dalmatian"` (wrapper to `r pkg("jags")` and `r pkg("nimble")`)
- **Model summary**: `r github("bbolker/broom.mixed")`, `r pkg("insight")`
- **Variable selection & Model Averaging**: `r pkg("LMERConvenienceFunctions")`, `r pkg("MuMIn")`
## Inference
### Hypothesis testing

- **fixed effects**: `r pkg("car")`, `r pkg("lmerTest")`, `r pkg("RVAideMemoire")`, `r pkg("emmeans")`, `r pkg("afex")`, `r pkg("pbkrtest")`, `r pkg("CLME")`, 
- **random effects**: `r pkg("varTestnlme")`, `r pkg("RLRsim")`, `r pkg("mvctm")`
### Prediction and Estimation

- `r pkg("emmeans")`, `r pkg("effects")`, `r pkg("margins")`, `r pkg("MarginalMediation")`
### Bootstrapping

- `r pkg("pbkrtest")`, `r pkg("lme4", priority = "core")` (`lme4::bootMer()` function), `r pkg("lmeresampler")`, `r pkg("glmmboot")`
### Power analysis

-`r pkg("longpower")`, `r pkg("clusterPower")`, `r pkg("powerlmm")`, `r pkg("pass.lme")`
## Other
### Commercial software interfaces

- [Mplus](https://www.statmodel.com/): `r pkg("MplusAutomation")` 
- [ASReml R](https://vsni.co.uk/software/asreml-r): `r pkg("asremlPlus")` 
- [Phoenix NLME software](http://www.certara.com/software/pkpd-modeling-and-simulation/phoenix-nlme/): `r pkg("Phxnlme")` 
### Links
#### Help

- [R-SIG-mixed-models mailing list](https://stat.ethz.ch/mailman/listinfo/r-sig-mixed-models) for discussion of mixed-model-related questions, course announcements, etc..
- [[r] + [mixed-models] tags on Stack Overflow](http://stackoverflow.com/questions/tagged/r+mixed-models)
- [Cross Validated](http://stats.stackexchange.com)
#### Other Software

- [ASReml-R](https://vsni.co.uk/software/asreml-r)
- [assist](https://yuedong.faculty.pstat.ucsb.edu/software.html)
- [INLA](http://www.r-inla.org/home)
- [JAGS](https://mcmc-jags.sourceforge.io/)
- [Stan](https://mc-stan.org) 
- [Zelig Project](https://zeligproject.org/)

#### Books

- *[SAS System for Mixed Models](https://v8doc.sas.com/sashtml/hrddoc/indfiles/55235.htm)*
- *[Generalized Linear Mixed Models](https://www.taylorfrancis.com/books/mono/10.1201/b13151/generalized-linear-mixed-models-walter-stroup)*
- *[Bayesian Data Analysis in Ecology using R, BUGS and Stan](https://www.elsevier.com/books/bayesian-data-analysis-in-ecology-using-linear-models-with-r-bugs-and-stan/korner-nievergelt/978-0-12-801370-0)* 
- *[Linear Mixed-Effects Models: A Step-by-Step Approach](https://link.springer.com/book/10.1007/978-1-4614-3900-4)* 
- *[Linear Mixed Models. An Introduction with applications in Veterinary Research](https://www.ilri.org/publications/linear-mixed-model-introduction-applications-veterinary-research)*

