

<!--    Submission: Open an issue in the ctv repository with the title CRAN task view proposal: MyTopic where MyTopic should be replaced with the intended name of the task view in CamelCase.
    Scope: Include one or two paragraphs about the scope of the task view, outlining the inclusion and exclusion criteria as well as relevant sections within the topic.
	-->
	
# Scope

*Mixed models* are a broad class of statistical models used to analyze data where observations can be assigned *a priori* to discrete groups, and where the parameters describing the differences between groups are treated as random variables. They are also described as *multilevel*, or *hierarchical*,  models; *longitudinal* data are often analyzed in this framework.  Mixed models can be fitted in either frequentist or Bayesian frameworks.

**Scope**: only including models that incorporate continuous (usually although not always Gaussian) latent variables; this excludes packages handling hidden Markov Models, finite (discrete) mixture models, latent Markov models, etc. We exclude general frameworks for implementing latent-variable models (Stan, TMB, greta, NIMBLE, etc.), but do include packages for fitting mixed models that are built on these platforms. 


<!--
    Packages: Include a tentative list of packages for the task view. This should encompass the "core" packages and a collection of relevant packages, ideally grouped by sections within the topic. It is not important, yet, that the list of sections or packages is already exhaustive.
    Overlap: Comment on potential overlap with already existing task views as well as with task views that might be created (or split off) in the future.
    Maintainers: The proposal should be made by the person willing to act as the principal maintainer for the task view. Furthermore, task views should have teams of additional 1-5 co-maintainers and possible candidates for these can be listed as well in the proposal.
-->
