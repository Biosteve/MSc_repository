gProfPlot2<-function(setofgenes, species, terms, top, header, plot=T, total=F){ 
  library(gprofiler2)
  library(ggplot2)
  # terms can be combinations of c("GO:BP", "GO:MF", "GO:CC", "KEGG", "HP", "HPA", "WP")
  if(terms=="GO"){terms1=c("GO:BP", "GO:MF", "GO:CC")}
  if(terms=="KEGG"){terms1="KEGG"}
  if(terms=="other"){terms1=c("HP", "HPA", "WP")}
  if(terms=="all"){terms1=c("GO:BP", "GO:MF", "GO:CC", "KEGG", "HP", "HPA", "WP")}
  gost(query=as.character(setofgenes), organism=species, significant=T, sources=terms1, evcodes=T)->out
  out<-out$result
  out<-out[which(out$term_size<=500 & out$term_size>=10),]
  out<-out[order(out$p_value/out$term_size, decreasing=F),]
  #out$new_term<-paste(out$source,":",out$term_name) # this version makes naming of terms more concise and short but may result in errors if many terms are returned
  out$new_term<-paste(out$term_id,":",out$term_name) # this version ensures no duplicates exist in the term factors
  out$new_term<-factor(out$new_term, levels = out$new_term[order(out$p_value, decreasing=T)])
    #
  if(plot==F & total==T) {return(out)}
  head(out,top)->out
  outlist<-data.frame("TermCategory"=out$source, "TermID"=out$term_id, "TermName"=out$term_name, "TermSize"=out$term_size, "log(p)"=-log10(out$p_value), "NumberOfGenes"=out$intersection_size, "Genes"=out$intersection)
  if(plot==F & total==F) {return(outlist)}
  #  
  funcplot<-ggplot(out, aes(x=new_term, y=-log10(p_value), label=round(-log10(p_value)),5)) + 
    geom_point(stat='identity', aes(col=-log10(p_value)), size=5)  +
    geom_text(color="black", size=2) +
    scale_color_gradient(low="steelblue4", high="firebrick4")+
    labs(title=header, 
         subtitle="Functional Enrichments") + 
    ylim(0, max(-log10(out$p_value))+1) +
    theme(text = element_text(size=10)) +
    coord_flip()
  if(plot==T){
  return(funcplot)
  }
}

# =============================================================================

gProfPlot2(as.character(gene_1304_list_sum$Gene), "hsapiens", top = 20,
           header="List 1304 genes (0C/ 1P/ 7C)", terms="all" ,
           plot = T,  total = F)

gProfPlot2(as.character(gene_2965_list_sum$Gene), "hsapiens", top = 20,
           header="List 2965 genes (0 c / 1 c / 7 c)", terms="all" ,
           plot = T,  total = F) 

gProfPlot2(as.character(gene_6193_list_sum$Gene), "hsapiens", top = 20,
           header="List 6193 genes (0 c / 1 c / 7 P)", terms="all" ,
           plot = T,  total = F) 

gProfPlot2(as.character(gene_11935_list_sum$Gene), "hsapiens", top = 20,
           header="List 11935 genes (0 C / 1 P / 7 P)", terms="all" ,
           plot = T,  total = F)

gProfPlot2(as.character(gene_26999_list_sum$Gene), "hsapiens", top = 20,
           header="List 26999 genes (0 P / 1 P / 7 P)", terms="all" ,
           plot = T,  total = F) 

