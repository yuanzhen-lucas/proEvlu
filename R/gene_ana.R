

#' Obtain the chromosomal and genetic information of the species from the protein file
#'
#' @param protein_df a data frame from function \code{\link{get_protein}} output
#' @param taxonomyID taxonomy ID of your species
#'
#' @return a list include two data frame(One contains chromosomal information
#'     and the other contains genetic information )
#'
#' @export
#'
#' @examples
#'
#'
#'
#' @author Yuan Zhen
#'
gene_sturcture = function(protein_df,taxonomyID){
  id_ref=entrez_link(dbfrom="taxonomy", id=taxonomyID, db="genome")

  genome_info=entrez_summary(db="genome",id=id_ref$links$taxonomy_genome)

  id_link_ref_gen=entrez_link(dbfrom = "genome",db="nucleotide",id=genome_info$uid)

  id_link_ref_gen=entrez_link(dbfrom = "genome",db="nucleotide",id=genome_info$uid)

  idd_test=list()
  if(length(id_link_ref_gen$links$genome_nuccore) < 10 ){

    web_load <- entrez_post(db="nucleotide", id=id_link_ref_gen$links$genome_nuccores)
    idd_nuccc=entrez_summary(db="nucleotide",web_history = web_load)

    idd_test=c(idd_test,idd_nuccc)

  } else {
    gene_caps = suppressWarnings(split(sample(id_link_ref_gen$links$genome_nuccore,replace = FALSE),1:(length(id_link_ref_gen$links$genome_nuccore) %/% 10)))
    for(i in seq(length(gene_caps))) {
      web_load <- entrez_post(db="nucleotide", id=gene_caps[[i]])
      idd_nuccc=entrez_summary(db="nucleotide",web_history = web_load)

      idd_test=c(idd_test,idd_nuccc)

    }
  }

  ttyz=lapply(1:length(idd_test),function(a){
    x=idd_test[[a]]
    if(x$sourcedb=="refseq" && x$genome=="chromosome"){
      chr_pos=which(str_detect(str_split(x$subtype,"[ | ]")[[1]],"chromosome"))
      chr_length=max(x$statistics$count[(which(x$statistics$type=="Length"))])
      data.frame(caption=x$caption,title=x$title,gi=x$gi,
                 length=chr_length,
                 chr=str_split(x$subname,"[|]")[[1]][chr_pos])
    }

  })

  ttyz=do.call(rbind,ttyz)


  id_gene=entrez_link(dbfrom="protein", id=protein_df$caption, db="gene")
  id_gene_mm=entrez_summary(db="gene",id=id_gene$links$protein_gene)

  id_vvv=entrez_link(dbfrom="gene", id=id_gene$links$protein_gene, db="nucleotide")

  id_vvv_sm=entrez_summary(db="nucleotide",id=id_vvv$links$gene_nuccore)



  ggyz=lapply(id_gene_mm,function(x){
    data.frame(name=x$name,chr=x$chromosome,
               start=x$genomicinfo$chrstart,
               end=x$genomicinfo$chrstop)
  })

  ggyz=do.call(rbind,ggyz)

  list(chromosome_info=ttyz,gene_info=ggyz)
}



#' Draw the position of genes on chromosomes
#'
#' @param gene_plot a list from function \code{\link{gene_sturcture}} output
#'
#' @return
#' @export
#'
#' @examples
#'
#'
plot_genes <- function(gene_plot) {
  ## Restrict to chromosomes that are in data
  ttyz=gene_plot$chromosome_info
  ggyz=gene_plot$gene_info
  chrs_in_data <-
    ttyz[ttyz$chr %in% ggyz$chr,]
  chr_order <- order(as.numeric(chrs_in_data$chr))

  ggplot() +
    geom_linerange(aes(x = chr,
                       ymin = 0,
                       ymax = length/1e6),
                   size = 8,
                   colour = heat.colors(length(chrs_in_data$chr)),
                   data = chrs_in_data) +
    geom_text_repel(aes(x = chr,
                        y = start/1e6,
                        label = name),
                    nudge_x = 0.33,
                    data = ggyz,
                    size=3,  fontface="bold",
                    box.padding=unit(0.5, "lines"), point.padding=unit(1.6, "lines"),
                    segment.color = "blue", segment.size = 0.5,
                    arrow = arrow(length=unit(0.01, "npc")),force = 1, max.iter = 3e3) +
    geom_point(aes(x = chr,
                   y = start/1e6,
                   label = name),data = ggyz,shape="-",size=4)+
    scale_y_reverse() +
    ## Fix ordering of chromosomes on x-axis
    scale_x_discrete(limits = chrs_in_data$chr[chr_order],
                     labels = paste0("chr",chrs_in_data$chr[chr_order])) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab("Chromosome") +
    ylab("Position (Mbp)")
  #annotate("text", x = chrs_in_data$chr, y = -3, label = paste0("chr",chrs_in_data$chr))

}
