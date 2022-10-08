

#' Get the pfm ID from the gene family name you gave
#'
#' @param protein_familyname The name of the gene family
#'
#' @return pfm IDs
#' @export
#'
#' @examples
#'
#'
#'
#' @seealso \code{\link{get_taxonomy}} which can be used to get taxonomy ID
#' @seealso \code{\link{get_protein}} which can be used to get taxonomy ID
#'
#' @author Zhen Yuan
#'
#'
get_pfm=function(protein_familyname){

  my_pfm=entrez_search(db="cdd", term =paste0(protein_familyname,"[all]"))
  my_sum=entrez_summary(db="cdd",id=my_pfm$ids)

  extract_name = as.character((extract_from_esummary(my_sum,c("uid","accession","title","subtitle","abstract","database"))))

  pfm_df = data.frame(matrix(unlist(extract_name),
                             nrow = length(extract_name)/6, byrow = T),
                      row.names = colnames(extract_name),
                      stringsAsFactors = F)

  colnames(pfm_df) = c("cdd_id","accession","title","subtitle","abstract","database")

  tax_id=entrez_link(dbfrom = "cdd",db="taxonomy",id=pfm_df$cdd_id)

  tax_summ=entrez_summary(db="taxonomy",id=tax_id$links$cdd_taxonomy)

  pfm_df$taxonomy_id=lapply(tax_summ,function(x)x$uid) %>% unlist
  pfm_df$division=lapply(tax_summ,function(x)x$division) %>% unlist
  pfm_df$scientificname=lapply(tax_summ,function(x)x$scientificname) %>% unlist
  pfm_df
}






#' Get the taxonomy ID from the species name you gave
#'
#' @param organism_name
#'
#' @return taxonomy IDs
#'
#'
#' @export
#'
#' @examples
#'
#' @seealso \code{\link{get_protein}} which can be used to get taxonomy ID
#' @seealso \code{\link{get_pfm}} which can be used to get pfm ID
#'
#' @author Zhen Yuan
#'
get_taxonomy=function(organism_name){
  my_taxonomy=entrez_search(db="taxonomy", term =paste0(organism_name,"[ORGN]"),retmax=60)
  my_taxonomy_id=entrez_summary(db="taxonomy",id=my_taxonomy$ids)


  extract_name = as.character((extract_from_esummary(my_taxonomy_id,c("uid","rank","division","scientificname"))))

  taxonomy_df = data.frame(matrix(unlist(extract_name),
                                  nrow = length(extract_name)/4, byrow = T),
                           row.names = colnames(extract_name),
                           stringsAsFactors = F)
  colnames(taxonomy_df)=c("taxonomy_id","rank","division","scientificname")
  taxonomy_df

}
