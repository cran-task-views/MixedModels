#!/usr/bin/env Rscript
# needs to adapted for our use cases

# Maintenance script to check CTV packages, URLs, and formatting.

library("ctv")
library("httr")
library("xml2")
library("magrittr")

ctvFile <- "MixedModels.md"
stopifnot(file.exists(ctvFile))

message("Building HTML and opening for viewing")
ctv::ctv2html(ctvFile)
htmlFile <- gsub(".md", ".html", ctvFile, fixed = TRUE)
browseURL(htmlFile)


message("Checking packages...")
packages <- check_ctv_packages(ctvFile)
packagesIssues <- lengths(packages) != 0
if (any(packagesIssues)) {
  warning("These packages need updating:", call. = FALSE, immediate. = TRUE)
  print(packages[packagesIssues])
}

message("Checking date...")
xml <- read_xml(htmlFile)
date_node <- xml_find_all(xml, "//meta[@name='DC.issued']")
cat(sprintf("Today is %s", Sys.Date()), "\n")
cat(sprintf("Task view last updated %s", xml_attr(date_node, "content")), "\n")
if (Sys.Date() != xml_attr(date_node, "content")) {
  warning("Don't forget to update the version", call. = FALSE, immediate. = TRUE)
}


message("Checking URLs...")

urls_all <- unique(xml_find_all(xml, "//a[@href]") %>% xml_attr(., "href")) 
urls <- urls_all[intersect(grep("^#.", urls_all, invert = TRUE),
                 grep("https://CRAN.R-project.org/.", urls_all, invert = TRUE))]
httr::set_config(timeout(1e6)) 
url_test <- rep(NA, length(urls))

## FIXME: progress bar? (back to a for loop?)
get_url <- function(url) {
    tt <- try(http_error(url,
                   config(ssl_verifypeer = 0L, ssl_verifyhost = 0L)),
              silent = TRUE)
    if (inherits(tt, "try-error")) return(NA)
    tt
}
url_test <- vapply(urls, get_url,
                   FUN.VALUE = logical(1))

                                        #url_test <- sapply(urls, try(http_error), config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))

## sometimes links come up error when they do work fine: false positive list
## (update as needed)

working_urls <- character(0)
bad_urls <- urls[url_test & !(urls %in% working_urls)]

if (length(bad_urls) > 0) {
    status <- vapply(bad_urls,
                     function(x) httr::GET(x)$status,
                     FUN.VALUE = integer(1))
    cat("Failed URLs:\n")
    vapply(status,
           function(x) http_status(x)$message,
           FUN.VALUE = character(1))
}

