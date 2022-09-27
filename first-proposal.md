

<!--    Submission: Open an issue in the ctv repository with the title CRAN task view proposal: MyTopic where MyTopic should be replaced with the intended name of the task view in CamelCase.
    Scope: Include one or two paragraphs about the scope of the task view, outlining the inclusion and exclusion criteria as well as relevant sections within the topic.
	-->
	
# Scope

*Mixed models* are a broad class of statistical models used to analyze data where observations can be assigned *a priori* to discrete groups, and where the parameters describing the differences between groups are treated as random variables. They are also described as *multilevel*, or *hierarchical*,  models; *longitudinal* data are often analyzed in this framework.  Mixed models can be fitted in either frequentist or Bayesian frameworks.

**Scope**: only including models that incorporate continuous (usually although not always Gaussian) latent variables; this excludes packages handling hidden Markov Models, finite (discrete) mixture models, latent Markov models, etc. We exclude general frameworks for implementing latent-variable models (Stan, TMB, greta, NIMBLE, etc.), but do include packages for fitting mixed models that are built on these platforms. 


## Packages

### Basic model fitting: (frequentist/LMM) `nlme`, `lme4`; (Bayesian/LMM) `MCMCglmm`, `rstanarm`, `brms`, `blme`; (frequentist/GLMM) `MASS`, `lme4`, `GLMMadaptive`, `hglm`; (Bayesian/GLMM) `MCMCglmm`, ...

## Overlap

There is some small overlap, e.g. with Bayesian Inference, DifferentialEquations, Spatial, but this is largely a distinct topic.

##  Maintainers

Julia Piaskowski

co-maintainer: Ben Bolker

