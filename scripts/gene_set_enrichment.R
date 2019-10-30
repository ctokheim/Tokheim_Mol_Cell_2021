suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(dplyr))
library(ReactomePA)
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(enrichplot))

# read data
znf2_df <- read.delim('/liulab/ctokheim/projects/degrader/transcription_factor/output/gata3_znf2/gene_correlations/BRCA.txt', sep='\t')
cterminal_df <- read.delim('/liulab/ctokheim/projects/degrader/transcription_factor/output/gata3_cterminal/gene_correlations/BRCA.txt', sep='\t')

# join two different mutations
df<-inner_join(znf2_df,cterminal_df, by="gene", suffix=c("_ZnF2", "_cterminal"))

# figure out what the entrez id is
tmp_gene_list <- df[,'gene']
geneToId <- bitr(tmp_gene_list, 
                 fromType="SYMBOL", 
                 toType="ENTREZID", 
                 OrgDb="org.Hs.eg.db")

# add in entrez id
df <- inner_join(df, geneToId, by=c('gene' = 'SYMBOL'))
df <- df[!is.na(df$ENTREZID),]

# read in msigdb genes
hallmark_path <- '/liulab/ctokheim/projects/degrader/transcription_factor/data/h.all.v6.2.entrez.gmt'
h <- read.gmt(hallmark_path)

#######################
# Run enrichment for cterminal
#######################

# prepare input
myGeneList <- df[,'GATA3_cterminal']
names(myGeneList) <- as.character(df[,'ENTREZID'])
myGeneList <- sort(myGeneList, decreasing = T)

# perform GSEA
hallmarkEnrich <- GSEA(myGeneList, TERM2GENE=h, verbose=FALSE, 
                       nPerm=100000, pvalueCutoff=0.1)

# convert to data frame
hallmark_enrich_df <- as.data.frame(hallmarkEnrich)

# save gsea plot
#gseaplot2(hallmarkEnrich, geneSetID = c('HALLMARK_ESTROGEN_RESPONSE_EARLY', 'HALLMARK_ESTROGEN_RESPONSE_LATE'),
          #pvalue_table = TRUE, color = c("#E495A5", "#7DB0DD"), subplots = 1)
#ggsave('gata3_cterminal.pdf')

# save gsea text results
write.table(hallmark_enrich_df, 'gata3_cterminal.txt', sep='\t', row.names=F)

#######################
# Run enrichment for ZnF2
#######################

# prepare input
myGeneList <- df[,'GATA3_ZnF2']
names(myGeneList) <- as.character(df[,'ENTREZID'])
myGeneList <- sort(myGeneList, decreasing = T)

# perform GSEA
hallmarkEnrich <- GSEA(myGeneList, TERM2GENE=h, verbose=FALSE, 
                       nPerm=100000, pvalueCutoff=0.1)

# convert to data frame
hallmark_enrich_df <- as.data.frame(hallmarkEnrich)

# save gsea plot
gseaplot2(hallmarkEnrich, geneSetID = c('HALLMARK_ESTROGEN_RESPONSE_EARLY', 'HALLMARK_ESTROGEN_RESPONSE_LATE'),
          pvalue_table = TRUE, color = c("#E495A5", "#7DB0DD"), subplots = 1)
ggsave('gata3_znf2.pdf')

# save gsea text results
write.table(hallmark_enrich_df, 'gata3_znf2.txt', sep='\t', row.names=F)


#######################
# Run enrichment for CTNNB1
#######################
# read data
ctnnb1_df <- read.delim('/liulab/ctokheim/projects/degrader/transcription_factor/output/gata3_znf2/gene_correlations/BRCA.txt', sep='\t')

# join two different mutations
df<-inner_join(znf2_df,cterminal_df, by="gene", suffix=c("_ZnF2", "_cterminal"))

# figure out what the entrez id is
tmp_gene_list <- df[,'gene']
geneToId <- bitr(tmp_gene_list, 
                 fromType="SYMBOL", 
                 toType="ENTREZID", 
                 OrgDb="org.Hs.eg.db")

# add in entrez id
df <- inner_join(df, geneToId, by=c('gene' = 'SYMBOL'))
df <- df[!is.na(df$ENTREZID),]


#######################
# Run enrichment for degron potential scores
#######################

degron_df <- read.delim('data/GPS/nterm_deepDegron/nterm_degron_predictions.txt')

# figure out what the entrez id is
tmp_gene_list <- degron_df[,'Gene_ID']
geneToId <- bitr(tmp_gene_list, 
                 fromType="SYMBOL", 
                 toType="ENTREZID", 
                 OrgDb="org.Hs.eg.db")

# add in entrez id
degron_df <- inner_join(degron_df, geneToId, by=c('Gene_ID' = 'SYMBOL'))
degron_df <- degron_df[!is.na(degron_df$ENTREZID),]

# prepare input
myGeneList <- degron_df[,'regulatory.potential']
names(myGeneList) <- as.character(degron_df[,'ENTREZID'])
myGeneList <- sort(myGeneList, decreasing = T)


keggEnrich <- gseKEGG(geneList     = myGeneList,
                     organism     = 'hsa',
                     nPerm        = 100000,
                     #minGSSize    = 10,
                     pvalueCutoff = 0.1,
                     verbose      = FALSE)
kegg_enrich_df <- as.data.frame(keggEnrich)

goEnrich2 <- gseGO(geneList     = myGeneList,
                   OrgDb        = org.Hs.eg.db,
                   keytype      = "ENTREZID",
                   ont          = "MF",
                   nPerm        = 100000,
                   pvalueCutoff = 0.1,
                   verbose      = FALSE
                   #maxGSSize = 1000
                   #minGSSize    = 100,
                  )
enrich_df <- as.data.frame(goEnrich2)
goEnrich3 <- simplify(goEnrich2, cutoff=0.7, by="p.adjust", select_fun=min)

reactomeEnrich <- gsePathway(geneList     = myGeneList,
                             organism     = 'human',
                             nPerm        = 100000,
                             #minGSSize    = 10,
                             pvalueCutoff = 0.1,
                             verbose      = FALSE)
reactome_enrich_df <- as.data.frame(reactomeEnrich)





require(magrittr) # because of the %<>% operator which is a pipe %>% with reassigning back `<-`
setMethod("simplify", signature(x="gseaResult"),
          function(x, cutoff=0.7, by="p.adjust", select_fun=min, measure="Wang", semData=NULL) {
            if (!x@setType %in% c("BP", "MF", "CC"))
              stop("simplify only applied to output from enrichGO...")
            x@result %<>% simplify_internal(., cutoff, by, select_fun,
                                            measure, x@setType, semData)
            return(x)
          }
)

# from:
# https://github.com/GuangchuangYu/clusterProfiler/blob/master/R/simplify.R
# I added `packagename::` in front of some function names, since this code is outside the package
# and does not "see" some of the packages/package functions imported 
# to the package environment
# but basically this is the unchanged `simplify_internal` function definition.

simplify_internal <- function(res, cutoff=0.7, by="p.adjust", select_fun=min, measure="Rel", ontology, semData) {
  if (missing(semData) || is.null(semData)) {
    if (measure == "Wang") {
      semData <- GOSemSim::godata(ont = ontology)
    } else {
      stop("godata should be provided for IC-based methods...")
    }
  } else {
    if (ontology != semData@ont) {
      msg <- paste("semData is for", semData@ont, "ontology, while enrichment result is for", ontology)
      stop(msg)
    }
  }
  
  sim <- GOSemSim::mgoSim(res$ID, res$ID,
                          semData = semData,
                          measure=measure,
                          combine=NULL)
  
  ## to satisfy codetools for calling gather
  go1 <- go2 <- similarity <- NULL
  
  sim.df <- as.data.frame(sim)
  sim.df$go1 <- row.names(sim.df)
  sim.df <- tidyr::gather(sim.df, go2, similarity, -go1)
  
  sim.df <- sim.df[!is.na(sim.df$similarity),]
  
  ## feature 'by' is attached to 'go1'
  sim.df <- merge(sim.df, res[, c("ID", by)], by.x="go1", by.y="ID")
  sim.df$go2 <- as.character(sim.df$go2)
  
  ID <- res$ID
  
  GO_to_remove <- character()
  for (i in seq_along(ID)) {
    ii <- which(sim.df$go2 == ID[i] & sim.df$similarity > cutoff)
    ## if length(ii) == 1, then go1 == go2
    if (length(ii) < 2)
      next
    
    sim_subset <- sim.df[ii,]
    
    jj <- which(sim_subset[, by] == select_fun(sim_subset[, by]))
    
    ## sim.df <- sim.df[-ii[-jj]]
    GO_to_remove <- c(GO_to_remove, sim_subset$go1[-jj]) %>% unique
  }
  
  res[!res$ID %in% GO_to_remove, ]
}

