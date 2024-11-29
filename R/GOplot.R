#' A simple function for Gene Ontology enrichment analysis and plotting.
#' @author qinti
#' @param gene
#' genelist from DEG
#' @return
#' ggplot2 plot object
#' @export
#' @examples
#' library(GOqt)
#' genelist = c("gene1","gene2","gene3")
#' GO_plot(genelist)+ggtitle("GO enrichment of diff_genes")
requireNamespace("org.Hs.eg.db", quietly = TRUE)
GO_plot = function(gene){
go_result <- clusterProfiler::enrichGO(gene = gene,
                   OrgDb=org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   minGSSize = 1,
                   pvalueCutoff = 0.99,
                   qvalueCutoff = 0.99,
                   readable = TRUE)
go_result_df =  go_result@result |>
    group_by(ONTOLOGY) |>
    arrange(p.adjust) |>
    slice_head(n=10) |>
    ungroup()
go_result_df = go_result_df[go_result_df$ONTOLOGY %in% c("BP","CC","MF"),]

mytheme <- theme(
axis.title = element_text(size = 13),
axis.text = element_text(size = 9),
plot.title = element_text(size = 14,
hjust = 0.5,
face = "bold"),
legend.title = element_text(size = 13),
legend.text = element_text(size = 11),
plot.margin = margin(t = 5.5,r = 10,l = 5.5,b = 5.5)
)
plot = ggplot(go_result_df, aes(Description, Count)) +
  geom_col(aes(fill = ONTOLOGY), width = 0.1) +
  geom_point(aes(size = -log10(p.adjust), color = ONTOLOGY)) +
  coord_flip() +
  labs(x = "", y = "Numbers of genes") +
  scale_color_manual(values = c(
    "#852f88",
    "#eb990c",
    "#0f8096"
  )) +
  scale_fill_manual(values = c(
    "#852f88",
    "#eb990c",
    "#0f8096"
  )) +
  scale_size_continuous(range = c(2, 5)) +
  theme_test()+mytheme+scale_x_discrete(labels = function(x) str_wrap(x, width = 70))+
  facet_grid(ONTOLOGY~.,scales = "free")+theme(strip.background = element_blank(),strip.text = element_blank())
return(plot)}


