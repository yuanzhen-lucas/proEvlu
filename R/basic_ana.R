



#' Protein families are obtained from NCBI by PFM number and taxonomy ID
#'
#' @param pfm pfm ID
#' @param taxonomyID taxonomy ID
#'
#' @return A data frame with protein accession ids ,descriptions,
#' sourcedb,taxonomy ID and organism
#' @export
#'
#' @examples
#'
#' \dontrun{
#'
#' my_pfm=get_pfm(protein_familyname="Dof")
#'
#' taxonomy_id=get_taxonomy(organism_name="soybean")
#'
#'###Please note that you may have multiple outputs for these two results, so choose the one that suits you best.
#'
#' protein_info=get_protein(pfm=my_pfm,taxonomyID=taxonomy_id)
#'
#' }
#'
#'
#' @seealso \code{\link{get_pfm}} which can be used to get pfm ID
#'
#' @seealso \code{\link{get_taxonomy}} which can be used to get taxonomy ID
#'
#'
#' @author Zhen Yuan
#'
get_protein=function(pfm="PF02701",taxonomyID="3983"){
  id_ref=entrez_link(dbfrom="taxonomy", id=taxonomyID, db="genome")
  if(is.null(id_ref$links$taxonomy_genome)){
    stop("Since the genome of the species you submitted could not be found in NCBI,
         we are sorry not to be able to perform a quick analysis for you")
  }

  id_cdd=entrez_search(db="cdd", term =paste0(pfm,"[all]"))
  id=entrez_link(dbfrom="cdd", id=id_cdd$ids, db="protein")
  protein_ids=id$links$cdd_protein
  extract_name=c()

  if(length(protein_ids) < 300 ){
    web_load <- entrez_post(db="protein", id=protein_ids)
    protein_summary = entrez_summary(db="protein", web_history=web_load)
    extract_name = as.character((extract_from_esummary(prot_summary,
                                                       c("caption","title",
                                                         "sourcedb", "taxid","organism"))))

  } else {
    protein_ids = suppressWarnings(split(sample(protein_ids,replace = FALSE),
                                         1:(length(protein_ids) %/% 280)))
    for(i in seq(length(protein_ids))) {
      web_load <- entrez_post(db="protein", id=protein_ids[[i]])
      protein_summary = entrez_summary(db="protein", web_history=web_load)
      extract_name = c(extract_name, as.character((extract_from_esummary(protein_summary,
                                                                         c("caption","title",
                                                                           "sourcedb", "taxid","organism")))))
    }
  }

  protein_df = data.frame(matrix(unlist(extract_name),
                                 nrow = length(extract_name)/5, byrow = T),
                          row.names = colnames(extract_name),
                          stringsAsFactors = F)

  colnames(protein_df) = c("caption","title","sourcedb", "taxid","organism")
  protein_df %>% filter(sourcedb=="refseq"&taxid==taxonomyID)

}
