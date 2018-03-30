# graph_plotter

graph_plotter <- function(overlap_matrix, allele, group_colors, suffix, cut,out_path){
  
  graph_network <- graph_from_adjacency_matrix(as.matrix(overlap_matrix), weighted = T, mode = "undirected")
  V(graph_network)
  if(length(E(graph_network))==0){
    print("No edges.")
    return(NA)
  }
  E(graph_network)
  # removing single parameters
  bad.vs<-V(graph_network)[graph_network$weight< 1]
  #graph_network_cut<-delete.vertices(graph_network, bad.vs)
  graph_network_cut <- graph_network
  
  
  #plot(graph_network_cut)
  pdf(file=paste0(out_path,"overlap_graph/",cut,"/", allele,"/",allele, suffix,".pdf"))
  par(mai=c(0,0,1,0)) 			#this specifies the size of the margins. the default settings leave too much free space on all sides (if no axes are printed)
  plot(graph_network_cut,				#the graph to be plotted
       layout=layout.fruchterman.reingold,	# the layout method. see the igraph documentation for details
       vertex.label.dist=0.25,			#puts the name labels slightly off the dots
       vertex.color = group_colors, #sapply(colnames(overlap_matrix_protein),function(x){ if(grepl("1",x)){return("green")}else if(grepl("2",x)){return("blue")}else{return("red")}}),
       vertex.frame.color=group_colors, 		#the color of the border of the dots 
       vertex.label.color='black',		#the color of the name labels
       vertex.label.font=0,			#the font of the name labels
       #vertex.label=V(bsk.network)$name,		#specifies the lables of the vertices. in this case the 'name' attribute is used
       vertex.label.cex=0.5,			#specifies the size of the font of the labels. can also be made to vary
       #vertex.shapes = 0.5
       vertex.size=3,
       edge.arrow.size = 0,
       edge.width=abs(E(graph_network_cut)$weight)*6
  )
  dev.off()
}