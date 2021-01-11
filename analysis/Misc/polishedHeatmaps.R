########### Ileum heatmap figures - Immune response genes #############
library(edgeR)
library(gplots)
library(RColorBrewer)
library(tidyverse)

## Load data
cnt1 <- read_delim("G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript/mallardLPAIV/extData/ileumCounts.txt", delim = "\t") %>% distinct(gene, .keep_all = TRUE)
covars <- read_delim("G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript/mallardLPAIV/extData/MALL70_SSgroup_RAW_11.5.20.csv", delim = ",")
annot <- read_delim("G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript/mallardLPAIV/extData/annotationGenes121320.txt", delim = "\t")


cnt <- cnt1 %>%
  as_tibble() %>%
  select(-gene)

newNames.ileum <- names(cnt) %>%
  as_tibble() %>%
  mutate(value = gsub("Iluem_", "", .$value)) %>%
  separate(value, into = c(NA, "pt1", "pt2")) %>%
  mutate(pt1 = na_if(pt1, "MALL"),
         pt2 = na_if(pt2, "ileum")) %>%
  mutate(bird = ifelse(is.na(pt1), pt2, pt1)) %>%
  mutate(bird = as.numeric(as.character(bird)),
         bird = formatC(bird, width = 2, format = "d", flag = "0")) %>%
  mutate(newName = paste0("ileum_", bird))

colnames(cnt) <- newNames.ileum$newName

cnt <- cnt[,order(colnames(cnt))] %>%
  mutate(gene = cnt1$gene) %>%
  column_to_rownames("gene") %>%
  as_tibble(rownames = NA)

### Covars
covars <- covars %>%
  mutate(bird = as.numeric(bird)) %>%
  mutate(bird = formatC(bird, width = 2, format = "d", flag = "0")) %>%
  arrange(bird) %>%
  mutate(bird = as.factor(bird)) %>%
  mutate(group = str_remove(group, "-"))

covars.ileum <- covars %>%
  filter(bird %in% newNames.ileum$bird)

cnt <- cnt %>% select(-ileum_01, -ileum_72)

covars.ileum <- covars.ileum %>% filter(bird != "01",
                                        bird != "72")

#Convert to DGEList object
dge <- DGEList(counts=cnt)
dge$genes <- annot[match(rownames(dge$counts), annot$ensembl_gene_id),]

#Assign groups in DGElist
group <- recode(covars.ileum$group, C1 = "Ctl", C29 = "Ctl") %>%
  as.factor()
dge$samples$group <- group

### Divide out I1 and I2
dge.I12 <- dge[,grep('\\bI1\\b|\\bI2\\b' , dge$samples$group)]
lcpm <- cpm(dge.I12, log = TRUE)




#### HEATMAP ####
targets = read_delim("../heatMapTargets.csv", delim = ",") %>%
  filter(tissue == "ileum")
# Get the gene names for DE genes
lcpm2 <- as_tibble(dge.I12$counts, rownames = NA) %>%
  cpm(., log = TRUE)

newNames <- colnames(lcpm2) %>%
  as_tibble() %>%
  separate(value, into = c(NA, "bird")) %>%
  left_join(covars) %>%
  unite("newNames", SSgroup.virus.sac, group, bird, remove = FALSE)

colOrder <- newNames %>%
  mutate(SSgroup.virus.sac = fct_relevel(SSgroup.virus.sac, "LOW", "MODERATE", "HIGH")) %>%
  arrange(group, SSgroup.virus.sac, bird) %>%
  unite("newNames", SSgroup.virus.sac, group, bird)

colnames(lcpm2) <- newNames$newNames
lcpm3 <- lcpm2[, colOrder$newNames]

heatmap.annot <- annot %>%
  filter(ensembl_gene_id %in% targets$ensembl_gene_id) %>%
  mutate(hgnc_symbol = replace_na(hgnc_symbol, ".")) %>%
  unite(annot.hm, c(ensembl_gene_id, hgnc_symbol), sep = ": ", remove = FALSE) %>%
  select(annot.hm, ensembl_gene_id) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

deGeneCounts <- lcpm3 %>%
  as_tibble() %>%
  mutate(ensembl_gene_id = rownames(dge.I12$counts)) %>%
  filter(ensembl_gene_id %in% targets$ensembl_gene_id) %>%
  left_join(heatmap.annot) %>%
  column_to_rownames("annot.hm") %>%
  select(-ensembl_gene_id) %>%
  as.matrix()

## Set up palette
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

hc <- hclust(as.dist(1-cor(t(deGeneCounts))))

# Plot the heatmap
heatmap.2(deGeneCounts,
          Colv = FALSE,
          Rowv = as.dendrogram(hc),
          col = rev(morecols(50)),
          trace = "none",
          colsep = c(5, 9, 12, 16, 20),
          dendrogram = "row",
          density.info = "none",
          key = FALSE,
          margins = c(10, 14),
          scale ="row")




############ BURSA #############
########### IImmune response genes #############
## Load data
cnt <- read.table("G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript/mallardLPAIV/extData/all.genes.results.csv", header = TRUE) %>% as_tibble(rownames = NA)
covars <- read_delim("G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript/mallardLPAIV/extData/MALL70_SSgroup_RAW_12.3.20.csv", delim = ",")
annot <- read_delim("G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript/mallardLPAIV/extData/annotationGenes121320.txt", delim = "\t")

newNames <- names(cnt) %>%
  as_tibble() %>%
  separate(value, into = c("pt1", "pt2", "pt3")) %>%
  select(pt2, pt3) %>%
  mutate(pt2 = as.numeric(pt2),
         pt2 = formatC(pt2, width = 2, format = "d", flag = "0")) %>%
  mutate(newName = paste0("bursa_", pt2))

colnames(cnt) <- newNames$newName

cnt <- cnt %>%
  select(order(colnames(cnt)))

covars <- covars %>%
  mutate(bird = as.numeric(bird)) %>%
  mutate(bird = formatC(bird, width = 2, format = "d", flag = "0")) %>%
  arrange(bird) %>%
  mutate(bird = as.factor(bird)) %>%
  mutate(group = str_remove(group, "-"))

## Remove birds with issues
cnt <- cnt %>% select(-bursa_01, -bursa_72)

covars <- covars %>% filter(bird != "01", bird != "72")

#Convert to DGEList object
dge <- DGEList(counts=cnt)
dge$genes <- annot[match(rownames(dge$counts), annot$ensembl_gene_id),]

#Assign groups in DGElist
group <- recode(covars$group, C1 = "Ctl", C29 = "Ctl") %>%
  as.factor()
dge$samples$group <- group

### Divide out I1 and I2
dge.I12 <- dge[,grep('\\bI1\\b|\\bI2\\b' , dge$samples$group)]
lcpm <- cpm(dge.I12, log = TRUE)




#### HEATMAP ####
targets = read_delim("../heatMapTargets.csv", delim = ",") %>%
  filter(tissue == "bursa")
# Get the gene names for DE genes
lcpm2 <- as_tibble(dge.I12$counts, rownames = NA) %>%
  cpm(., log = TRUE)

newNames <- colnames(lcpm2) %>%
  as_tibble() %>%
  separate(value, into = c(NA, "bird")) %>%
  left_join(covars) %>%
  unite("newNames", SSgroup.virus.sac, group, bird, remove = FALSE)

colOrder <- newNames %>%
  mutate(SSgroup.virus.sac = fct_relevel(SSgroup.virus.sac, "LOW", "MODERATE", "HIGH")) %>%
  arrange(group, SSgroup.virus.sac, bird) %>%
  unite("newNames", SSgroup.virus.sac, group, bird)

colnames(lcpm2) <- newNames$newNames
lcpm3 <- lcpm2[, colOrder$newNames]

heatmap.annot <- annot %>%
  filter(ensembl_gene_id %in% targets$ensembl_gene_id) %>%
  mutate(hgnc_symbol = replace_na(hgnc_symbol, ".")) %>%
  unite(annot.hm, c(ensembl_gene_id, hgnc_symbol), sep = ": ", remove = FALSE) %>%
  select(annot.hm, ensembl_gene_id) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE)

deGeneCounts <- lcpm3 %>%
  as_tibble() %>%
  mutate(ensembl_gene_id = rownames(dge.I12$counts)) %>%
  filter(ensembl_gene_id %in% targets$ensembl_gene_id) %>%
  left_join(heatmap.annot) %>%
  column_to_rownames("annot.hm") %>%
  select(-ensembl_gene_id) %>%
  as.matrix()

## Set up palette
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

hc <- hclust(as.dist(1-cor(t(deGeneCounts))))

# Plot the heatmap
heatmap.2(deGeneCounts,
          Colv = FALSE,
          Rowv = as.dendrogram(hc),
          col = rev(morecols(50)),
          trace = "none",
          colsep = c(6, 12, 15, 19, 23),
          dendrogram = "row",
          density.info = "none",
          key = FALSE,
          margins = c(10, 24),
          scale ="row")
