library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(ggplot2)
library(dplyr)

#https://yulab-smu.top/biomedical-knowledge-mining-book/clusterProfiler-dplyr.html

setwd("./collaboration/ARIC_Cancer/all/original/pathway")
load("./collaboration/ARIC_Cancer/all/original/results/cox_fit_fullmodel_liver_anno.RData")
entrez_ids = res_final$`Entrez Gene ID`[which(res_final$fdr < 0.05)]
uniprot_id = res_final$uniprot_id[which(res_final$fdr < 0.05)]

# -------- KEGG Pathway Enrichment --------
kegg_enrich <- enrichMKEGG(
  gene = entrez_ids,
  #universe      = res_final$uniprot_id,
  #keyType = 'uniprot',
  organism = 'hsa',
  pvalueCutoff = 1
)


go_enrich <- enrichGO(gene   = entrez_ids,
        # universe      = names(geneList),
         OrgDb         = org.Hs.eg.db,
        # keyType = 'uniprot',
         ont           = "CC",
         pAdjustMethod = "BH",
         pvalueCutoff  = 1,
         qvalueCutoff  = 0.05,
         readable      = TRUE)



go_enrich_bp <- enrichGO(gene   = entrez_ids,
                      # universe      = names(geneList),
                      OrgDb         = org.Hs.eg.db,
                      # keyType = 'uniprot',
                      ont           = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 1,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)




go_enrich_mf <- enrichGO(gene   = entrez_ids,
                         # universe      = names(geneList),
                         OrgDb         = org.Hs.eg.db,
                         # keyType = 'uniprot',
                         ont           = "MF",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 1,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)


wp <- enrichWP(entrez_ids, organism = "Homo sapiens") 

library(ReactomePA)
reactome <- enrichPathway(gene=entrez_ids, pvalueCutoff = 0.05, readable=TRUE)


save(reactome,go_enrich,go_enrich_bp,go_enrich_mf,wp,kegg_enrich,file="pathway_liver.RData")

# Prepare data for plots
#kegg_df <- as.data.frame(kegg_enrich)
load("./collaboration/ARIC_Cancer/all/original/pathway/pathway_liver.RData")
go_df <- as.data.frame(go_enrich)
go_enrich_bp_df <- as.data.frame(go_enrich_bp)
go_enrich_mf_df <- as.data.frame(go_enrich_mf)
#wp_df <- as.data.frame(wp)
reactome_df <-as.data.frame(reactome)

reactome_df$pathway = "Reactome"
go_enrich_bp_df$pathway = "GO Biological Process"
go_enrich_mf_df$pathway = "GO Molecular Function"
go_df$pathway = "GO Cellular Component"

full_results_print = rbind(reactome_df,go_enrich_bp_df,go_enrich_mf_df,go_df)
full_results_print$GeneRatio = as.character(full_results_print$GeneRatio)
full_results_print$GeneRatio <- ifelse(grepl("^\\d+/\\d+$", full_results_print$GeneRatio),
                          paste0("'", full_results_print$GeneRatio),
                          full_results_print$GeneRatio)

#full_results_print$GeneRatio <- paste0("='", full_results_print$GeneRatio, "'")
write.csv(full_results_print,file="table_print.csv")

