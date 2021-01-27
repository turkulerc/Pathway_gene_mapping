
geneNumberInGO <- function(genename, GOname){
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    install.packages("org.Hs.eg.db")
  if (!requireNamespace("AnnotationDbi", quietly = TRUE))
    install.packages("AnnotationDbi")
  if (!requireNamespace("GO.db", quietly = TRUE))
    install.packages("GO.db")
  if (!requireNamespace("tidyverse", quietly = TRUE))
    install.packages("tidyverse")
  require(org.Hs.eg.db)
  require(AnnotationDbi)
  require(GO.db)
  require(tidyverse)
  obj <- "org.Hs.egGO2ALLEGS"
  orgPkg <- "org.Hs.eg.db"
  GeneID.PathID <- toTable(tryCatch(getFromNamespace(obj,orgPkg), error=function(e) FALSE))
  annots_gene <- select(org.Hs.eg.db, keys=GeneID.PathID$gene_id,columns=c("SYMBOL","GENENAME"), keytype="ENTREZID") %>% 
    dplyr::filter(SYMBOL %in% genename) %>% 
    dplyr::distinct(SYMBOL, ENTREZID, .keep_all = TRUE)
  annots_GO <- Term(GeneID.PathID$go_id) 
  annots_GO <- tibble::tibble(go_id = names(annots_GO),
                   GO_name = annots_GO) %>% 
    dplyr::distinct(go_id, GO_name, .keep_all = TRUE) %>% 
    dplyr::filter(GO_name %in% GOname)
  tmp <- dplyr::left_join(annots_gene, GeneID.PathID, by = c("ENTREZID" = "gene_id")) %>% 
    dplyr::left_join(annots_GO) %>% 
    dplyr::select(SYMBOL, ENTREZID, go_id, GO_name) %>% 
    dplyr::distinct(SYMBOL, GO_name, .keep_all = TRUE) %>%
    tidyr::drop_na() %>% 
    dplyr::group_by(GO_name) %>% 
    dplyr::summarise(Ngenes = n())
}
  
geneNumberInKEGG <- function(genename, KEGGname){
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    install.packages("org.Hs.eg.db")
  if (!requireNamespace("AnnotationDbi", quietly = TRUE))
    install.packages("AnnotationDbi")
  if (!requireNamespace("tidyverse", quietly = TRUE))
    install.packages("tidyverse")
  if (!requireNamespace("RCurl", quietly = TRUE))
    install.packages("RCurl")
  require(org.Hs.eg.db)
  require(AnnotationDbi)
  require(tidyverse)
  require(RCurl)
  annots_KEGG <- limma::getGeneKEGGLinks("hsa", convert=FALSE)
  annots_gene <- select(org.Hs.eg.db, keys=GeneID.PathID$gene_id,columns=c("SYMBOL","GENENAME"), keytype="ENTREZID") %>% 
    dplyr::filter(SYMBOL %in% genename) %>% 
    dplyr::distinct(SYMBOL, ENTREZID, .keep_all = TRUE)
  uniqueKEGG <- tibble::tibble(KEGGid = unique(annots_KEGG$PathwayID)) %>% 
    dplyr::mutate(KEGGid = str_replace(KEGGid, "path:", ""), KEGGname = NA)
  for(i in 1:nrow(uniqueKEGG)){
    uniqueKEGG$KEGGname[i] <- RCurl::getURL(paste0("http://togows.dbcls.jp/entry/pathway/", uniqueKEGG$KEGGid[i], "/name"))
  }
  
  
  
  require(rvest)
  for(length(unique(annots_KEGG$PathwayID))){
    kegg <- read_html(paste0("https://www.kegg.jp/dbget-bin/www_bget?pathway+", annots_KEGG$PathwayID), encoding = "ISO-8859-1")  
    tmp <- kegg %>% html_nodes("td") %>% html_text()