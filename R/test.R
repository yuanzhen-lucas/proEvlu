
library(rentrez)

library(tidyverse)

library(Biostrings)

library(ggrepel)



prot=test(id$links$cdd_protein_summary)
prot=read_rds("prot.Rds")

prot_f=prot %>% filter(sourcedb=="refseq"&organism=="pepper")

id_prot=entrez_summary(db="protein",id=prot_f$caption)
all_prot <- entrez_fetch(db="protein", id=prot_f$caption, rettype="fasta")

prot_seq=str_split(all_prot,"\n")[[1]] 
prot_seq=prot_seq[!str_count(prot_seq)==0]
str_name=which(str_detect(prot_seq,">"))
str_q=lapply(seq_along(str_name),function(x){
  if(x==length(str_name)){
    str_c(prot_seq[(str_name[x]+1):(length(prot_seq))],collapse = "")
  }else{
    str_c(prot_seq[(str_name[x]+1):(str_name[x+1]-1)],collapse = "")
  }
  
})%>%unlist


names(str_q)=prot_seq[str_name]%>%str_remove_all(">")

pt_seq=AAStringSet(str_q)






id_gene=entrez_link(dbfrom="protein", id=prot_f$caption, db="gene")
id_gene_mm=entrez_summary(db="gene",id=id_gene$links$protein_gene)







dna_seq_ids <- entrez_link(dbfrom="gene", id=id_gene$links$protein_gene, db="nuccore")
yz <- entrez_summary(db="gene",id=id_gene$links$protein_gene)

all_seq <- entrez_fetch(db="gene", id=id_gene$links$protein_gene, rettype="fasta")

write_rds(all_seq,"dna_seq.Rds")


dna_seq=str2_sting(all_seq)





########## plot gene structure #######









#####protein##
library(protr)
library(Biostrings)
library(msa)
library(seqinr)
library(tinytex)
library(ape)

id_pro_mm=entrez_summary(db="protein",id=prot_f$caption)

pro_seq=entrez_fetch(db="protein",id=prot_f$caption[1:10],rettype="fasta")
pro_seq2=entrez_fetch(db="protein",id=prot_f$caption[11:24],rettype="fasta")


pro1=str2_sting(pro_seq,type="pro")
pro2=str2_sting(pro_seq2,type="pro")
pro=c(pro1,pro2)









#######  get protein ID from NCBI #####


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
    protein_ids = suppressWarnings(split(sample(protein_ids,replace = FALSE),1:(length(protein_ids) %/% 280)))
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


test=get_protein()


#######  get gene sturcture ######

gene_sturcture = function(protein_df=test,taxonomyID="3983"){
  id_ref=entrez_link(dbfrom="taxonomy", id=taxonomyID, db="genome")
  
  genome_info=entrez_summary(db="genome",id=id_ref$links$taxonomy_genome)
  
  id_link_ref_gen=entrez_link(dbfrom = "genome",db="nucleotide",id=genome_info$uid)
  
  idd_nuccc=entrez_summary(db="nucleotide",id=id_link_ref_gen$links$genome_nuccore)
  
  ttyz=lapply(1:length(idd_nuccc),function(a){
  
    x=idd_nuccc[[a]]
    if(x$sourcedb=="refseq" && x$genome=="chromosome"){
      chr_pos=which(str_detect(str_split(x$subtype,"[ | ]")[[1]],"chromosome"))
      chr_length=max(x$statistics$count[(which(x$statistics$type=="Length"))])
      data.frame(caption=x$caption,title=x$title,gi=x$gi,
                 length=chr_length,
                 chr=str_split(x$subname,"[ | ]")[[1]][chr_pos])
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


############ plot prootein evlution  ####

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




protein_properties <- function(pro_seq){
  data.frame(aminoacid_length=lapply(pro_seq,function(x)Peptides::lengthpep(x)) %>%unlist(),
             isoelectic_point=lapply(pro_seq,function(x)Peptides::pI(x)) %>%unlist(),
             molecular_weight_kDa=lapply(pro_seq,function(x)Peptides::mw(x)) %>%unlist() / 1000,
             theoretical_net_charge=lapply(pro_seq,function(x)Peptides::charge(x)) %>%unlist(),
             hydrophobicity_index=lapply(pro_seq,function(x)Peptides::hydrophobicity(x)) %>%unlist(),
             instability_index=lapply(pro_seq,function(x)Peptides::instaIndex(x)) %>%unlist()
             )
}


pro_tree = function(pro_seq,plot=FALSE){
  
  pro_ali=msa(pro_seq)
  sdist <- stringDist(as(pro_ali, "BStringSet"), method="hamming")
  clust <- hclust(sdist, method = "single")
 
  if(!plot){
    plot(clust)
  }

  
}




yz = pro_seq %>% as.character()








