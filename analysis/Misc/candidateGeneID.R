##### Query annotations to identify candidate genes ######
library(tidyverse)
candidateList <- read_delim("G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript/CandidateGeneList.csv", delim = ",")

annot <- read_delim("G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript/mallardLPAIV/extData/annotationTrans103120.txt", delim = "\t")

for(term in candidateList$term) {
  assign(paste0(term, ".df"),
         annot %>%
           filter_all(any_vars(str_detect(., term))) %>%
         mutate(searchTerm = term))
}

rm(candidateList, annot)
DF_obj <- lapply(ls(), get)

candidates <- Reduce(rbind, DF_obj) %>%
  distinct(ensembl_transcript_id, .keep_all = TRUE)

write.csv(candidates, "G:/Shared drives/RNA Seq Supershedder Project/MALL DE manuscript/candidatesMallardTranscripts.csv")
