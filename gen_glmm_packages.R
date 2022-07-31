## utilities/misc code for finding interesting packages related to mixed models
library(tidyverse) ## general 
library(miniCRAN)
library(crandep)
library(igraph)
library(ctv)
## library(packdep) ## archived ...
library(packageRank)

options(repos = c(CRAN = "https://cloud.r-project.org"))

## 1. plot reverse dependencies of lme4
## dependency types to include
## A depends on B: B automatically gets loaded
## A imports B: A uses functions from B, B *must* be installed
## A suggests B: A optionally uses functions from B
## A enhances B: the authors of the package think these go together
## "reverse-depends": reverse dependencies of lme4 = all the packages that depend on lme4
rd <- c("Reverse depends", "Reverse imports", "Reverse linking to", "Reverse suggests")
dd <- get_dep("lme4",rd)
## plot(igraph::graph_from_data_frame(dd))  ## too many!
## nrow(dd) ## 425, 31 July 2022
## grep/regular expression tricks

## 2. find many (not all) GLMM packages
a1 <- available.packages()
## grep("lmm",rownames(a1),value=TRUE,ignore.case=TRUE)
## lm followed by ("m or e") followed by (a character not "t" or the end of the string)
regexps <- c("lm(m|e([^t]|$))")  ## was using "mixed" but ... ? should check,
find_pkgs <- function(x) grep(x, rownames(a1), value=TRUE, ignore.case=TRUE)
regex_pkgs <- character(0)
for (r in regexps) {
    regex_pkgs <- union(regex_pkgs, find_pkgs(r))
}
## false pos
fpos <- c("palmerpenguins", "curtailment", "yamlme", "mailmerge")


## false negatives: some known-interesting pkgs *not* picked up by regex
## (check MixedModels.ctv for some more)
## fneg <- c("SASmixed", "broom.mixed",
##           "pbkrtest", "emmeans", "mgcv", "gamm4",
##           "brms", "rstanarm", "pez", "merDeriv", "repeated", "hglm",
##           "geesmv", "geepack", "influence.ME", "cAIC4", "HLMdiag", "lmmfit",
##           "iccbeta", "DHARMa", "effects", "rockchalk",
##           "arm", "performance", "car",
##           "ez", "afex", "RVAideMemoire", "geoRglm", "GLMMarp", "spaMM",
##           "polytomous", "ordinal", "longpower")


regex_pkgs <- setdiff(regex_pkgs, fpos)
rr <- read.ctv("MixedModels.md")
bad_pkgs <- unlist(check_ctv_packages("MixedModels.md"))
ctv_pkgs <- setdiff(rr$packagelist[,"name"], bad_pkgs)

omit_pkgs <- "MASS"
## don't want to include MASS in ranking (it's only in there for glmmPQL)

focal_pkgs <- union(ctv_pkgs, regex_pkgs) |> setdiff(omit_pkgs)

length(focal_pkgs) ## 128

## now extract
pkg_rd <- (expand_grid(name=focal_pkgs, type=rd))

## clunky
pb <- txtProgressBar(max=nrow(pkg_rd), style=3)
i <- 0
ff <- function(name,type) {
    ## cat(".")
    ## cat(name, "\n")
    i <<- i+1
    setTxtProgressBar(pb,i)
    gd <- get_dep(name,type)
    res <- tibble(focal=name, type, from = gd$from, to = gd$to)
    if (nrow(res) == 0) return(NULL)
    return(res)
}


if (file.exists("all_deps.rds")) {
    all_deps <- readRDS("all_deps.rds")
} else {
    ## THIS BIT IS SLOW, WATCH OUT ... (~ 6 minutes)
    system.time(all_deps <- (pkg_rd
        |> pmap(ff)
        |> bind_rows()
    )
    )
    saveRDS(all_deps, "all_deps.rds")
}

## disregard numbers, for purposes of plotting
unique_deps <- (all_deps
    |> select(focal, from, to)
    |> unique()
    |> filter(to %in% focal_pkgs)
)

rdg1 <-igraph::graph_from_data_frame(unique_deps[,c("to", "from")])

## plot(rdg1)
     
## 3. collect importance measures

cc <- eigen_centrality(rdg1)$vector
central_tbl <- tibble(focal=names(cc), central=cc)

a1a <- (all_deps
    |> mutate_at("type", str_remove, "Reverse ")
    |> mutate_at("type",
                  ~ case_when(. %in% c("depends", "imports") ~ "strong",
                              TRUE ~ "weak"))
)

a1 <- (a1a
    |> drop_na()
    |> count(focal, type)
)

a_tot <- (a1
    |> group_by(focal)
    |> summarise(n=sum(n),.groups="drop")
    |> mutate(type="total")
)

a2 <- (bind_rows(a1,a_tot)
    |> pivot_wider(names_from=type, values_from=n, values_fill=0)
    ## restore packages with no depends
    |> full_join(tibble(focal=focal_pkgs), by="focal")
    |> mutate(across(where(is.integer), replace_na, 0L))
)

## 
## SLOW the first time
pp <- packageRank(focal_pkgs)
pp2 <- (pp$package.data
    |> as_tibble()
    |> select(focal=packages,downloads,percentile)
)

a3 <- (full_join(pp2, a2, by="focal")
    |> full_join(central_tbl, by="focal")
    |> mutate_at("central",replace_na,0)
    |> mutate(score=percentile/100+strong/max(strong)+central, .before = downloads)
    |> mutate(across(where(is.double), round, 3))
    |> mutate(in_taskview = focal %in% ctv_pkgs, .before = downloads)
    |> rename(package = "focal")
    |> arrange(desc(score))
)

## a3
## View(a3)
## add descriptions??

write_csv(a3, "glmm_packages.csv")
rmarkdown::render("glmm_packages_meta.rmd", output_format = "md_document")
