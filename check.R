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
  warning("The packages need updated", call. = FALSE, immediate. = TRUE)
  print(packages[packagesIssues])
}


message("Checking URLs...")

xml <- read_xml(htmlFile)
urls_all <- unique(xml_find_all(xml, "//a[@href]") %>% xml_attr(., "href")) 
urls <- urls_all[intersect(grep("^#.", urls_all, invert = TRUE),
                 grep("https://CRAN.R-project.org/.", urls_all, invert = TRUE))]
httr::set_config(timeout(1e6)) 
url_test <- rep(NA, length(urls))
for(i in 1:length(urls)) {
  url_test[i] = try(http_error(urls[i], config(ssl_verifypeer = 0L, ssl_verifyhost = 0L)))
}

url_test <- sapply(urls, try(http_error), config(ssl_verifypeer = 0L, ssl_verifyhost = 0L))

# sometimes links come up error when they do work fine: 
  # (update as needed)
working_urls <- NA

cat("here are the URL's that did failed the test:\n")
urls[url_test == "TRUE" & !(urls %in% working_urls)] # %>% datapasta::vector_paste()
cat("here are the URL's with error messages:\n")
urls[grep("^Error.", url_test)] 


message("Checking Date...")
date_node <- xml_find_all(xml, "//meta[@name='DC.issued']")
sprintf("Today is %s", Sys.Date())
sprintf("This was last updated %s", xml_attr(date_node, "content"))
if (Sys.Date() != xml_attr(date_node, "content")) {
  warning("Don't forget to update the version", call. = FALSE, immediate. = TRUE)
}


