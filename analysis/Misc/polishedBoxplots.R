########### Publication-ready boxplots ############
targets = c("ENSAPLG00020015991",
            "ENSAPLG00020017110")
## Libraries
library(edgeR)
library(gridExtra)
library(tidyverse)

### Get gene level counts
## Load and prep
cnt.gene <- read.table("G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript/mallardLPAIV/extData/all.genes.results.csv", header = TRUE) %>% as_tibble(rownames = NA)
annot.gene <- read_delim("G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript/mallardLPAIV/extData/annotationGenes103120.txt", delim = "\t")
cnt.trans <- read.table("G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript/mallardLPAIV/extData/all.isoforms.results.csv", header = TRUE) %>% as_tibble(rownames = NA)
annot.trans <- read_delim("G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript/mallardLPAIV/extData/annotationTrans103120.txt", delim = "\t")
covars <- read_delim("G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript/mallardLPAIV/extData/MALL70_SSgroup_RAW_11.5.20.csv", delim = ",")

newNames <- names(cnt.gene) %>%
  as_tibble() %>%
  separate(value, into = c("pt1", "pt2", "pt3")) %>%
  select(pt2, pt3) %>%
  mutate(pt2 = as.numeric(pt2),
         pt2 = formatC(pt2, width = 2, format = "d", flag = "0")) %>%
  mutate(newName = paste0("bursa_", pt2))

colnames(cnt.gene) <- newNames$newName
colnames(cnt.trans) <- newNames$newName


cnt.gene <- cnt.gene %>%
  select(order(colnames(cnt.gene)))
cnt.trans <- cnt.trans %>%
  select(order(colnames(cnt.trans)))


covars <- covars %>%
  mutate(bird = as.numeric(bird)) %>%
  mutate(bird = formatC(bird, width = 2, format = "d", flag = "0")) %>%
  arrange(bird) %>%
  mutate(bird = as.factor(bird)) %>%
  mutate(group = str_remove(group, "-")) %>%
  filter(bird != "01", bird != "72")

## Remove birds with issues
cnt.gene <- cnt.gene %>% select(-bursa_01, -bursa_72)
cnt.trans <- cnt.trans %>% select(-bursa_01, -bursa_72)

#Convert to DGEList object
dge.gene <- DGEList(counts=cnt.gene)
dge.trans <- DGEList(counts=cnt.trans)

#CPM and log-CPM
lcpm.gene <- cpm(dge.gene, log = TRUE) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column("id")

lcpm.trans <- cpm(dge.trans, log = TRUE) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column("id")



plot.tib <- lcpm.trans %>%
  bind_rows(lcpm.gene) %>%
  pivot_longer(cols = starts_with("bursa"), names_to = "sample", values_to = "lcpm") %>%
  filter(id %in% targets) %>%
  separate(sample, into = c("tissue", "bird")) %>%
  mutate(bird = as.factor(bird)) %>%
  left_join(covars, by = "bird") %>%
  select(id, bird, lcpm, group) %>%
  mutate(group = recode(group, C29 = "Control"),
         group = recode(group, C1 = "Control"),
         tmpID = ifelse(group == "Control", 'Control', 'Infected'),
         id = recode(id, ENSAPLT00020026836 = "HSP8-Trans",
                     ENSAPLG00020017110 = "HSPA8-Gene/Trans",
                     ENSAPLG00020015991 = "FOSB-Gene"))


plot.tib %>%
  mutate(group = fct_relevel(group, "Control", "I1", "I2", "I5", "I15", "I29")) %>%
  ggplot(aes(x = group, y = lcpm, fill = id)) +
  facet_grid(. ~ tmpID, scales = "free", space = "free") +
  ylab("Log2(Counts per million)") +
  xlab("Group") +
  geom_boxplot() +
  theme_classic(base_size = 18) +
  scale_fill_grey(
    start = 0.5,
    end = 0.95,
    na.value = "red",
    aesthetics = "fill"
  ) +
  theme(legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank())
