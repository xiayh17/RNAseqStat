# enhance_enrichGO <- function(...,gene,
#                              orgDb = 'org.Hs.eg.db',
#                              ont="BP",
#                              keyType = "SYMBOL",
#                              pvalueCutoff=0.01,
#                              qvalueCutoff=0.2) {
#   ego <- enrichGO(gene = gene,
#            OrgDb = orgDb,
#            keyType = keyType,
#            ont = ont,
#            pvalueCutoff = pvalueCutoff,
#            pAdjustMethod = "BH",
#            qvalueCutoff = qvalueCutoff,
#            minGSSize = 10,
#            maxGSSize = 500,
#            readable = FALSE,
#            pool = FALSE,...)
#
#   return(ego)
# }
#
# enhance_enrichGO()
#
# data(geneList, package = "DOSE")
# de <- names(geneList)[1:100]
# DEG_df_g <- cut_much(DEG_df,x = "log2FoldChange",y = "pvalue",cut_FC = 2,cut_P = 0.01)
# gene_down <- row.names(DEG_df_g[which(DEG_df_g$group == "Down"),])
# yy <- enhance_enrichGO(gene = gene_down, orgDb = 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01)
#
# enrichplot::cnetplot(yy)
# head(yy)




