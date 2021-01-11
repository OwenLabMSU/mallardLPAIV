library(tidyverse)
dat <- read_delim("./Read_counts.filtered.txt",
                  delim = "\t",
                  col_names = FALSE)
dat <- dat %>%
  mutate(variable = rep(c("sample", "readCount"), nrow(dat) / 2),
         key = rep(1:(nrow(dat) / 2), each = 2)) %>%
  pivot_wider(id_cols = key, names_from = variable, values_from = X1) %>%
  select(-key) %>%
  filter(!grepl("excess", sample)) %>%
  mutate(readCount = as.numeric(readCount))


dat1 <- dat %>%
  filter(!grepl("MALL", sample)) %>%
  separate(sample, into = c("bird", "tissue", NA, NA, NA, NA))


dat2 <- dat %>%
  filter(grepl("MALL", sample)) %>%
  separate(sample, into = c(NA, "tissue", "bird", NA, NA, NA)) %>%
  mutate(tissue = "ileum")

left_join(dat1, dat2) %>%
  group_by(tissue) %>%
  summarize(minReads = min(readCount),
            maxReads = max(readCount),
            meanReads = mean(readCount))
