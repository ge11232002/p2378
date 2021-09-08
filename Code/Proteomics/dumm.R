library(tidyverse)

exdata <- readRDS(url("https://github.com/tidyverse/tidyr/files/3918551/exdata.zip"))
res <- tidyr::complete(
  exdata$data,
  tidyr::nesting(!!!syms(exdata$set1)),
  tidyr::nesting(!!!syms(exdata$set2))
)

exdata <- readRDS(url("https://github.com/tidyverse/tidyr/files/3918551/exdata.zip"))
res <- tidyr::complete(
  exdata$data,
  tidyr::nesting(!!!syms(exdata$set1)),
  tidyr::nesting(!!!syms(exdata$set2[c(5,1:4)]))
)
