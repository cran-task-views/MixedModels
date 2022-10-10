## FIXME: set up proper make rule for all_deps.rds
glmm_packages:
	R CMD BATCH --vanilla gen_glmm_packages.R

taskview:
	Rscript -e "library(ctv); ctv2html(read.ctv('MixedModels.md'))"
