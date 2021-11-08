# 24hrs -------------------------------------------------------------------
library(tidyverse)
dataset <- read_tsv("/Volumes/Transcend/p2378/dataset.tsv") %>%
  filter(str_detect(Name, "24hrs")) %>%
  mutate("R1" = str_replace(`Read1 [File]`, fixed("p2378/NovaSeq_20211020_NOV982_o25922_DataDelivery/"), "")) %>%
  select(Name, R1)

md5sum <- read_delim("/Volumes/Transcend/p2378/md5sums.txt", col_names = FALSE)
ans <- dataset %>% left_join(md5sum, by = c("R1"="X3"))
write_tsv(ans, "GEO.txt")

library(readxl)
counts <- read_excel("/Volumes/Transcend/p2378_24h/Count_QC-raw-count.xlsx")
counts <- counts %>% select(`Feature ID`, contains("24hrs"))
write_tsv(counts, "Count_QC-raw-count.txt")


# 6hrs --------------------------------------------------------------------
library(tidyverse)
dataset <- read_tsv("/Volumes/Transcend/p2378/dataset.tsv") %>%
  filter(str_detect(Name, "6hrs")) %>%
  mutate("R1" = str_replace(`Read1 [File]`, fixed("p2378/NovaSeq_20211020_NOV982_o25922_DataDelivery/"), "")) %>%
  select(Name, R1)

md5sum <- read_delim("/Volumes/Transcend/p2378/md5sums.txt", col_names = FALSE)
ans <- dataset %>% left_join(md5sum, by = c("R1"="X3"))
write_tsv(ans, "GEO.txt")

md5sum <- read_lines("/Volumes/Transcend/p2378/md5sums.txt")
write_lines(str_subset(md5sum, "6hrs"), "/Volumes/Transcend/p2378/md5sums_6hrs.txt")

library(readxl)
counts <- read_excel("/Volumes/Transcend/p2378/Count_QC-raw-count.xlsx")
counts <- counts %>% select(`Feature ID`, contains("6hrs"))
write_tsv(counts, "Count_QC-raw-count.txt")
