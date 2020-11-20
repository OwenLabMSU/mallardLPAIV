library(rtracklayer)

setwd("G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript")

cnt.gene <- read.table("./mallardLPAIV/extData/all.genes.results.csv", header = TRUE)
cnt.trans <- read.table("./mallardLPAIV/extData/all.isoforms.results.csv", header = TRUE)

ensIDs.gene <- rownames(cnt.gene) %>%
  as_tibble() %>%
  separate(value, into = c("ensID", NA, NA))

ensIDs.trans <- rownames(cnt.trans) %>%
  as_tibble() %>%
  separate(value, into = c("ensID", NA, NA))

library(biomaRt)
ensembl = useMart("ensembl",
                  dataset="aplatyrhynchos_gene_ensembl")

annotGene <- getBM(attributes=c('ensembl_gene_id', 'gene_biotype', 'hgnc_symbol', 'description', 'entrezgene_id'),
                   filters = 'ensembl_gene_id',
                   values = ensIDs.gene$ensID,
                   mart = ensembl)

annotTrans <- getBM(attributes=c('ensembl_transcript_id', 'transcript_biotype', 'hgnc_symbol', 'description', 'entrezgene_id'),
                   filters = 'ensembl_transcript_id',
                   values = ensIDs.trans$ensID,
                   mart = ensembl)


write.table(annotGene, file = "annotationGenes103120.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(annotTrans, file = "annotationTrans103120.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

'%ni%' <- Negate('%in%')
ensIDs %>% filter(ensID %ni% annotGene$ensembl_gene_id)




#### Ileum, reconstruct ensembl IDs
# Retrieve ensembl IDs for genes annotated with just hgnc symbol
# Use list of ensembl IDs to query the mart and make annotation file

library(biomaRt)
library(tidyverse)

cntTmp <- read_delim("G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript/adjusted_read_counts.csv", delim = ',') %>%
  as_tibble(rownames = NA) %>%
  select(gene) %>%
  mutate(ensembl_gene_id = ifelse(str_detect(gene, "ENS"), gene, NA)) %>%
  mutate(hgnc_symbol = ifelse(str_detect(gene, "ENS"), NA, gene))

hgnc.target <- cntTmp %>% drop_na(hgnc_symbol)


ensembl = useMart("ensembl",
                  dataset="aplatyrhynchos_gene_ensembl")

annotGene1 <- getBM(attributes=c('ensembl_gene_id', 'gene_biotype', 'hgnc_symbol', 'description', 'entrezgene_id'),
                    filters = 'hgnc_symbol',
                    values = hgnc.target$hgnc_symbol,
                    mart = ensembl) %>%
  select(hgnc_symbol, ensembl_gene_id)

cnt <- read_delim("G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript/adjusted_read_counts.csv", delim = ',') %>%
  as_tibble(rownames = NA) %>%
  select(-chr) %>%
  mutate(ensembl_gene_id = ifelse(str_detect(gene, "ENS"), gene, NA)) %>%
  mutate(hgnc_symbol = ifelse(str_detect(gene, "ENS"), NA, gene))


cnt.ens <- left_join(cnt, annotGene1, by = "hgnc_symbol") %>%
  rename(ensembl_gene_id = ensembl_gene_id.y) %>%
  mutate(gene = ifelse(str_detect(gene, "ENS"), gene, ensembl_gene_id)) %>%
  select(-ensembl_gene_id.x, -hgnc_symbol, -ensembl_gene_id) %>%
  drop_na(gene)


annotGene.ileum <- getBM(attributes=c('ensembl_gene_id', 'gene_biotype', 'hgnc_symbol', 'description', 'entrezgene_id'),
                     filters = 'ensembl_gene_id',
                     values = cnt.ens$gene,
                     mart = ensembl)

write.table(annotGene.ileum, file = "annotationGenesIleum110120.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(cnt.ens, file = "ileumCounts.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

