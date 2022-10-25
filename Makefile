## FIXME: set up proper make rule for all_deps.rds
glmm_packages:
	R CMD BATCH --vanilla gen_glmm_packages.R

taskview:
	Rscript -e "ctv::ctv2html(ctv::read.ctv('MixedModels.md'))"

check:
	Rscript check.R
	Rscript -e "ctv::check_ctv_packages('MixedModels.md')"
