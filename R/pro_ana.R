
#' Turn NCBI obtained information into DNAString or AAString objects
#'
#' @param seq
#' @param type
#'
#' @return
#' @export
#'
#' @examples
str2_sting=function(seq,type="dna"){
  prot_seq=str_split(seq,"\n")[[1]]
  prot_seq=prot_seq[!str_count(prot_seq)==0]
  str_name=which(str_detect(prot_seq,">"))
  str_q=lapply(seq_along(str_name),function(x){
    if(x==length(str_name)){
      str_c(prot_seq[(str_name[x]+1):(length(prot_seq))],collapse = "")
    }else{
      str_c(prot_seq[(str_name[x]+1):(str_name[x+1]-1)],collapse = "")
    }

  })%>%unlist


  names(str_q)=prot_seq[str_name]%>%str_remove_all(">") %>%
    str_split(" ") %>%
    map(.,function(x)x[[1]]) %>% unlist

  if(type=="dna"){
    DNAStringSet(str_q)
  }else{

  }
  AAStringSet(str_q)

}



#' Title
#'
#' @param protein_df
#'
#' @return
#' @export
#'
#' @examples
protein_seq <- function(protein_df=test){
  protein_cap=protein_df$caption
  pro_aa_seq=Biostrings::AAStringSet()
  if(length(protein_cap) < 10 ){
    web_load <- entrez_post(db="protein", id=protein_cap)

    pro_aa_seq=entrez_fetch(db="protein",web_history=web_load,rettype="fasta")


  } else {
    protein_caps = suppressWarnings(split(sample(protein_cap,replace = FALSE),1:(length(protein_cap) %/% 10)))
    for(i in seq(length(protein_caps))) {
      web_load <- entrez_post(db="protein", id=protein_caps[[i]])

      pro_seq=entrez_fetch(db="protein",web_history=web_load,rettype="fasta")

      pro_aa_seq=c(pro_aa_seq,str2_sting(pro_seq,type="pro"))

    }
  }
  pro_aa_seq

}


pro_seq=protein_seq()




#' Title
#'
#' @param pro_seq
#'
#' @return
#' @export
#'
#' @examples
protein_properties <- function(pro_seq){
  data.frame(aminoacid_length=lapply(pro_seq,function(x)Peptides::lengthpep(x)) %>%unlist(),
             isoelectic_point=lapply(pro_seq,function(x)Peptides::pI(x)) %>%unlist(),
             molecular_weight_kDa=lapply(pro_seq,function(x)Peptides::mw(x)) %>%unlist() / 1000,
             theoretical_net_charge=lapply(pro_seq,function(x)Peptides::charge(x)) %>%unlist(),
             hydrophobicity_index=lapply(pro_seq,function(x)Peptides::hydrophobicity(x)) %>%unlist(),
             instability_index=lapply(pro_seq,function(x)Peptides::instaIndex(x)) %>%unlist()
  )
}


#' Title
#'
#' @param pro_seq
#' @param plot
#'
#' @return
#' @export
#'
#' @examples
pro_tree = function(pro_seq,plot=FALSE){

  pro_ali=msa(pro_seq)
  sdist <- stringDist(as(pro_ali, "BStringSet"), method="hamming")
  clust <- hclust(sdist, method = "single")

  if(!plot){
    plot(clust)
  }


}
