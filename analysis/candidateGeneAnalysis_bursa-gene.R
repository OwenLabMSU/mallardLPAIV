## Load packages and data
library(edgeR)
library(nlme)
library(broom)
library(tidyverse)

setwd("G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript/mallardLPAIV/")

## Load data
annot <- read.delim("./extData/annotationGenes103120.txt", header = TRUE, sep = "\t")
cnt <- read.table("./extData/all.genes.results.csv", header = TRUE)
covars <- read.csv("./extData/MALL70_SSgroup_RAW_11.5.20.csv", header = TRUE)
targets <- read_delim("./extData/MALL_candidateGeneList_filtered.csv", delim = ",") %>%
  filter(include == "yes")


##### Prepare data ####
newNames <- names(cnt) %>%
  as_tibble() %>%
  separate(value, into = c("pt1", "pt2", "pt3")) %>%
  select(pt2, pt3) %>%
  mutate(pt2 = as.numeric(pt2),
         pt2 = formatC(pt2, width = 2, format = "d", flag = "0")) %>%
  mutate(newName = paste0("bursa_", pt2))

colnames(cnt) <- newNames$newName

cnt <- cnt %>%
  select(order(colnames(cnt))) %>%
  select(-bursa_01, -bursa_72)

### covars
covars <- covars %>%
  mutate(bird = as.numeric(bird)) %>%
  mutate(bird = formatC(bird, width = 2, format = "d", flag = "0")) %>%
  arrange(bird) %>%
  mutate(bird = as.factor(bird)) %>%
  mutate(group = str_remove(group, "-")) %>%
  filter(bird != "01", bird != "72")

#### Calculate log(CPM) and assemble master DFs ####
# Convert to DGEList object
dge <- DGEList(counts=cnt)
lcpm <- cpm(dge, log = TRUE)


# Make lcpm tibble for analysis
lcpm.all <- lcpm %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column("gene") %>%
  filter(gene %in% targets$ensembl_gene_id) %>%
  pivot_longer(cols = contains("_"),
               names_to = "sample",
               values_to = "lcpm") %>%
  separate(sample, into = c("tissue", "bird")) %>%
  mutate(bird = factor(bird)) %>%
  left_join(covars, by = "bird") %>%
  mutate(identifier = gene) %>%
  select(identifier, bird, tissue, lcpm, virus.sac, SSgroup.virus.sac, group, age, sex, Bursa_pool, Ileum_pool) %>%
  mutate(pool = ifelse(tissue == "bursa", Bursa_pool, Ileum_pool)) %>%
  mutate(log.virus.sac = log(virus.sac+0.01)) %>%
  mutate(group = recode(group,
                        C1 = "Ctl",
                        C29 = "Ctl")) %>%
  mutate(group = factor(group,
                        levels = c("Ctl", "I1", "I2", "I5", "I15", "I29")))

#### Gene query database ####
annot <- annot %>%
  mutate(hgnc_symbol = ifelse(hgnc_symbol == "", ".", hgnc_symbol))


##### Analysis function ######
candidateGeneAnalysis.Reg <- function(target, ...) {
  library(tidyverse)
  library(nlme)

  tryCatch({
    filtered.dat <- lcpm.all %>%
      filter(identifier == target)

    tryCatch({
      sum.overall <- filtered.dat %>%
        lme(lcpm ~ log.virus.sac + age + factor(sex), random = ~1|pool, data = ., control = lmeControl(opt = "optim")) %>%
        anova(.)},
      error=function(err) NA)


    ifelse(exists("sum.overall"),
           overall.result <- as_tibble(sum.overall$`p-value`[2]) %>%
             mutate(identifier = target,
                    subset = "Overall"),
           overall.result <- data.frame("value" = NA, "identifier" = NA, "subset" = NA))


    tryCatch({
      sum.I1 <- filtered.dat %>%
        filter(group == "I1") %>%
        lme(lcpm ~ log.virus.sac + age + factor(sex), random = ~1|pool, data = ., control = lmeControl(opt = "optim")) %>%
        anova(.)},
      error=function(err) NA)


    ifelse(exists("sum.I1"),
           I1.result <- as_tibble(sum.I1$`p-value`[2]) %>%
             mutate(identifier = target,
                    subset = "I1"),
           I1.result <- data.frame("value" = NA, "identifier" = NA, "subset" = NA))

    tryCatch({
      sum.I2 <- filtered.dat %>%
        filter(group == "I2") %>%
        lme(lcpm ~ log.virus.sac + age + factor(sex), random = ~1|pool, data = ., control = lmeControl(opt = "optim")) %>%
        anova(.)},
      error=function(err) NA)


    ifelse(exists("sum.I2"),
           I2.result <- as_tibble(sum.I2$`p-value`[2]) %>%
             mutate(identifier = target,
                    subset = "I2"),
           I2.result <- data.frame("value" = NA, "identifier" = NA, "subset" = NA))


    tryCatch({
      sum.I5 <- filtered.dat %>%
        filter(group == "I5") %>%
        lme(lcpm ~ log.virus.sac + age + factor(sex), random = ~1|pool, data = ., control = lmeControl(opt = "optim")) %>%
        anova(.)},
      error=function(err) NA)


    ifelse(exists("sum.I5"),
           I5.result <- as_tibble(sum.I5$`p-value`[2]) %>%
             mutate(identifier = target,
                    subset = "I5"),
           I5.result <- data.frame("value" = NA, "identifier" = NA, "subset" = NA))


    tryCatch({
      sum.I15 <- filtered.dat %>%
        filter(group == "I15") %>%
        lme(lcpm ~ log.virus.sac + age + factor(sex), random = ~1|pool, data = ., control = lmeControl(opt = "optim")) %>%
        anova(.)},
      error=function(err) NA)


    ifelse(exists("sum.I15"),
           I15.result <- as_tibble(sum.I15$`p-value`[2]) %>%
             mutate(identifier = target,
                    subset = "I15"),
           I15.result <- data.frame("value" = NA, "identifier" = NA, "subset" = NA))

    bind_rows(overall.result, I1.result, I2.result, I5.result, I15.result)
  }, error=function(e){})
}


#### Run analysis loop ####
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl)

set <- lcpm.all %>%
  group_by(identifier) %>%
  summarize(varLCPM = round(var(lcpm), 5)) %>%
  filter(varLCPM > 0)

finalMatrix <- foreach(z = unique(set$identifier), .combine = rbind) %dopar% {
  tmpMatrix = candidateGeneAnalysis.Reg(z)
  tmpMatrix
}

#stop cluster
stopCluster(cl)



## Identify significant patterns using the overall dataset
finalMatrix.sig <- finalMatrix %>%
  filter(subset == "Overall") %>%
  mutate(adj.p.value = p.adjust(value, method='fdr', n = nrow(.))) %>%
  filter(adj.p.value < 0.05) %>%
  mutate(comparison = "BG")


### Clean up and save ###
remove(
  annot,
  finalMatrix.BG,
  finalMatrix.BT,
  finalMatrix.IG,
  finalMatrix.IT,
  finalMatrix.BG.sig,
  finalMatrix.BT.sig,
  finalMatrix.IG.sig,
  finalMatrix.IT.sig,
  cnt.bursa.gene,
  cnt.bursa.trans,
  cnt.gene,
  cnt.ileum.gene,
  cnt.ileum.trans,
  cnt.trans,
  dge.bursa.gene,
  dge.bursa.trans,
  dge.ileum.gene,
  dge.ileum.trans,
  lcpm.bursa.gene,
  lcpm.bursa.trans,
  lcpm.bursa.tmp,
  lcpm.ileum.gene,
  lcpm.ileum.tmp,
  lcpm.ileum.trans,
  lcpm.trans,
  set,
  targets
)

save.image("candidateGeneAnalysis_regression-Ctl.Rws")
