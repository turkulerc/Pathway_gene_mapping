
# Provide GO term names in the example format - 
# spaces as word delimiters and without "GO": "B cell activation"
geneNumberInGO <- function(genename, GOname){
  if (is.null(genename)) stop("There is not any gene name that provided as input")
  if (is.null(GOname)) stop("There is not any GO term name that provided as input")
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
  GOname <- tolower(GOname)
  GeneID.PathID <- toTable(tryCatch(getFromNamespace(obj,orgPkg), error=function(e) FALSE))
  annots_gene <- select(org.Hs.eg.db, keys=GeneID.PathID$gene_id,columns=c("SYMBOL","GENENAME"), keytype="ENTREZID") %>% 
    dplyr::filter(SYMBOL %in% genename) %>% 
    dplyr::distinct(SYMBOL, ENTREZID, .keep_all = TRUE)
  annots_GO <- Term(GeneID.PathID$go_id) 
  annots_GO <- tibble::tibble(go_id = names(annots_GO),
                   GO_name = annots_GO) %>% 
    dplyr::distinct(go_id, GO_name, .keep_all = TRUE) %>% 
    dplyr::mutate(GO_name = tolower(GO_name)) %>% 
    dplyr::filter(GO_name %in% GOname)
  tmp <- dplyr::left_join(annots_gene, GeneID.PathID, by = c("ENTREZID" = "gene_id")) %>% 
    dplyr::left_join(annots_GO) %>% 
    dplyr::select(SYMBOL, ENTREZID, go_id, GO_name) %>% 
    dplyr::distinct(SYMBOL, GO_name, .keep_all = TRUE) %>%
    tidyr::drop_na() 
  table_tmp <- dplyr::group_by(tmp, GO_name) %>% 
    dplyr::mutate(genes = paste0(SYMBOL, collapse = ", ")) %>% 
    dplyr::select(GO_name, genes) %>% 
    dplyr::distinct(GO_name, genes, .keep_all = TRUE)
  table <-
    dplyr::group_by(tmp, GO_name) %>% 
    dplyr::summarise(Ngenes = n())
  return(list(genes = table_tmp, numbers = table))
}
  
# Provide KEGG pathway names in the example format - 
# all in small letters and spaces as word delimiters: "citrate cycle (tca cycle)"
geneNumberInKEGG <- function(genename, KEGGname){
  if (is.null(genename)) stop("There is not any gene name that provided as input")
  if (is.null(KEGGname)) stop("There is not any KEGG pathway name that provided as input")
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
  KEGGname <- tolower(KEGGname)
  annots_KEGG <- limma::getGeneKEGGLinks("hsa", convert=FALSE) %>% 
    dplyr::mutate(PathwayID = str_replace(PathwayID, "path:", ""))
  annots_gene <- select(org.Hs.eg.db, keys=annots_KEGG$GeneID,columns=c("SYMBOL","GENENAME"), keytype="ENTREZID") %>% 
    dplyr::filter(SYMBOL %in% genename) %>% 
    dplyr::distinct(SYMBOL, ENTREZID, .keep_all = TRUE)
  uniqueKEGG <- tibble::tibble(KEGGid = unique(annots_KEGG$PathwayID)) %>% 
    dplyr::mutate(KEGGid = str_replace(KEGGid, "path:", ""), KEGG_name = NA)
  for(i in 1:nrow(uniqueKEGG)){
    uniqueKEGG$KEGGname[i] <- RCurl::getURL(paste0("http://togows.dbcls.jp/entry/pathway/", uniqueKEGG$KEGGid[i], "/name"))
  }
  uniqueKEGG <- dplyr::mutate(uniqueKEGG, KEGG_name = str_replace(KEGG_name, " - Homo sapiens \\(human\\)\n", "")) %>% 
    dplyr::mutate(KEGG_name = tolower(KEGG_name))
  tmp <- dplyr::left_join(annots_gene, annots_KEGG, by = c("ENTREZID" = "GeneID")) %>% 
    dplyr::left_join(uniqueKEGG, by = c("PathwayID" = "KEGGid")) %>% 
    dplyr::select(SYMBOL, ENTREZID, PathwayID, KEGG_name) %>% 
    dplyr::distinct(SYMBOL, KEGG_name, .keep_all = TRUE) %>%
    tidyr::drop_na() %>% 
    dplyr::filter(SYMBOL %in% genename, KEGG_name %in% KEGGname)
  table_tmp <- dplyr::group_by(tmp, KEGG_name) %>% 
    dplyr::mutate(genes = paste0(SYMBOL, collapse = ", ")) %>% 
    dplyr::select(KEGG_name, genes) %>% 
    dplyr::distinct(KEGG_name, genes, .keep_all = TRUE)
  table <-
    dplyr::group_by(tmp, KEGG_name) %>% 
    dplyr::summarise(Ngenes = n())
  return(list(genes = table_tmp, numbers = table))
}