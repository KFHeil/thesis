#the following script allows to read edge lists and protein tables with additional vertex information

#remove all objects
rm(list=ls())

#we import the igraph library
library(igraph)
#library needed to fill matrix
library(plyr)

#############################################################

#we generate the plots based on the spectral algorithm

#function to set the graph layout
network_layout <- function(algrth, network, name_out){
  #we set structural properties of the graph, per subgraph (and later on plot the networks)
  if(algrth == "spectral"){
    for(b in 1:max(vertex_attr(network)$community)){
      #print(c("i", b))
      #print(c("alg$algorithm", alg$algorithm))
      sg <- induced_subgraph(network, V(network)[V(network)$community==b])
      lsg <- layout_(sg, in_circle()) * length(V(network)$community[V(network)$community == b])/1.5
      #modify the layout of the graph based on the lcf - kamada kawai computation
      # lg[V(network)$community == b,] <- t(t(lsg)+lcg[b+1,])
      lg[V(network)$community == b,] <- t(t(lsg)+lcg[b,])
    }
  } else {
    print("no")
    for(i in 1:length(algrth)){
      #we check if the algorithm is the spectral one - if so this functions changes slightly 
      # print(c("i", i))
      #    print(c("alg$algorithm", alg$algorithm))
      sg <- induced_subgraph(network, V(network)[algrth$membership==i])
      lsg <- layout_(sg, in_circle()) * sizes(algrth)[i]/1.5
      #modify the layout of the graph based on the lcf - kamada kawai computation
      # lg[algrth$membership==i,] <- t(t(lsg)+lcg[i+1,])
      #not entirely sure but needed change... also above
      lg[algrth$membership==i,] <- t(t(lsg)+lcg[i,])
    }
  }
  ########
  #SETTING GRAPH ATTRIBUTES
  #we set the graph properties (both of the below approach are working)
  #we set the size of the vertexex/nodes
  vert_lab_dist <- 0.3
  vert_lab_font <- 3
  
  #set the general color and size for the vertex labels
  V(network)$vertex.label.color="grey"
  V(network)$vertex.label.cex=0.5
  
  ########
  #computing the eigenvector centrality score to determint the size of the vertices
  eigen_centrality <- eigen_centrality(network)
  #the values are stored in
  #print(c(eigen_centrality$vector, "eigen_centrality$vector"))
  #we multiply the eigenvalue centrality values before shaping the network vertices
  eigen_centrality_mod <- eigen_centrality$vector * 10
  #and set the attribute
  V(network)$vertex.size = eigen_centrality_mod
  
  #plot the network to pdf
  plot(network, vertex.label = V(network)$gene_name, edge.color = E(network)$edge.color, layout = lg, vertex.label.color = V(network)$label.color, vertex.size=5, vertex.label.cex = 0.7, main=paste(folder, "g_biggest_component_", name_out, sep = ""), vertex.shape= V(network)$vertex.shapes, vertex.color = V(network)$color, sub=paste("vertices: ", gorder(network), " edges: ", nrow(as_edgelist(network)), " PD-hits: ", sum(vertex_attr(network)$PD), sep = ""))
  plot(network, main=paste(folder, "graph-", name_out, sep = ""), layout = lg, vertex.label.dist = vert_lab_dist, vertex.label.font = vert_lab_font, vertex.label.color = V(network)$vertex.label.color, vertex.size = V(network)$vertex.size, vertex.label.cex =V(network)$vertex.label.cex)
  write(paste("Number of PD-associated genes in the biggest connected component is: ", sum(vertex_attr(network)$PD), sep = " "), file=log_file, append = TRUE)
  #  dev.off()
}


#this functions generates output files with all the gene aliases in one row that are localised in one community
dataframe_communities <- function(network, name_test){
  #we generate output files with the proteins per community of each clustering algorithm
  #we generate an empty matrix
  #rows = number of communities
  #columns = max number of proteins per community
  rows = unique(vertex_attr(network)$community)
  #we determine the max number of proteins/nodes in the biggest cluster
  #we iter over the different clusters and in case the cluster size is bigger we update the value of ncols
  ncols = 0
  for(a in unique(vertex_attr(network)$community)){
    nodes_community = length(vertex_attr(network)$community[vertex_attr(network)$community == a])
    #    print(c(a, nodes_community))
    if(nodes_community < ncols){
      next
    } else {
      ncols = nodes_community
    }
  }
  #we generate an overview table
  obj=matrix(nrow=length(rows), ncol=ncols+1)
  #we iter over the number of communities per network
  for(t in 1:length(rows)){
    #  print(c("i", i))
    #  print(groups(alg)[[i]])
    #and add the protein identifiers to the matrix
    #we do this by selecting all the proteins in the community of interest
    obj[t,] <-c(t(rbind.fill.matrix(rows[t], vertex_attr(g_biggest_component)$name[vertex_attr(network)$community == rows[t]], rep(NA, ncols - length((vertex_attr(network))[vertex_attr(network)$community == rows[t]])))))
    #we add the data to the dataframe
  }
  #we create the folder for the community output files
  file_path = paste(homedir, "/", folder, sep="")
  #we finally print the data to an output file
  output_file = paste(file_path, name_test, "/communities_members.csv", sep = "")
  #print(c("output_file", output_file))
  #we set the row names to false - can be reset - but keeps community numbers of the whole graph before selecting the biggest component (not per default numbered from 1 - x, but can have random numbers)
  write.table(obj, output_file, sep = ',', row.names = F, col.names=F, na="", quote = FALSE)
}


#function to generate dataframe with proteins of interest in communities --> this can be used to run the hypergeometric enrichment test
dataframe_proteins <- function(gbc, name_test){
  #generate dataframe to fill with information of interest
  dataf <- data.frame(matrix(nrow=length(vertex_attr(gbc)$name), ncol=3))
  names(dataf)<-c("GeneID", "PD", "community")
  #we fill the data frame
  dataf$GeneID <- vertex_attr(gbc)$name
  dataf$community <- vertex_attr(gbc)$community
  dataf$PD <- vertex_attr(gbc)$PD
  print(dataf)
  #we generate the path to the output folder
  file_path_df = paste(homedir, "/", folder, sep="")
  write.csv(dataf, file = paste(file_path_df, name_test, "/data_frame.csv", sep=""), row.names = FALSE, quote = FALSE)
  cat(paste("data frame generated:", file_path_df, name_test, "/data_frame.csv", sep=""))
  community_statistics(dataf, name_test)
}

community_statistics <- function(dataframe, name_test){
  #we generate a matrix with one row per community in the network 
  obj <- data.frame(matrix(nrow=length(unique(dataframe$community)), ncol = 4))
  #print(obj)
  #possibly add PD_overlap
  names(obj)<-c("Community", "community_size", "PD", "PD_per_cent")
  #iter over the number of communities and calculate the respective occurences of PD hits
  #we additionally save the names of the hits in lists
  #we generate a dataframe to add associated proteins
  #we set the number of columns to the maximal of PD associated nodes in the network 
  obj_G <- matrix(nrow=length(unique(dataframe$community)), ncol=sum(dataframe$PD == 1)+1)
  names(obj_G) <- c("community", "PD")
  #we add a counter for rows in the dataframe
  counter = 0
  for(r in unique(dataframe$community)){
    #we add one to the counter for every row in the dataframe/community in the network
    counter = counter + 1
    #    print(c("counter", counter))
    #    print(c("r", r))
    PD <- 0
    PD_hits <- c()
    total <- 0
    for(u in 1:nrow(dataframe)){
      #      print(c("u", u))
      if(dataframe$community[u] == r){
        #        print(c("dataframe$community[u]", dataframe$community[u]))
        #        print(c("r", r))
        total <- total + 1
        #        print(c("community ", r))
        if(dataframe$PD[u] == 1){
          PD <- PD + 1
          PD_hits <- c(PD_hits, dataframe$GeneID[u])
        }
      }
    }
    obj$Community[counter] <- r
    obj$community_size[counter] <- total
    obj$PD[counter] <- PD
    obj$PD_per_cent[counter] <- (PD/total)*100
    
    #we add the PD associated nodes to the dataframes of interest
    #we do this by selecting all the proteins in the community of interest
    #we add the community number - associated nodes/hits and NAs in the remaining columns
    #    print("fill specific objects")
    obj_G[counter,] <-c(t(rbind.fill.matrix(r, as.numeric(PD_hits), rep(NA, sum(dataframe$PD == 1) - length(PD_hits)))))
    #    print(c("obj_G", obj_G))
  }
  #we generate the path to the output folder
  file_path_df = paste(homedir, "/", folder, sep="")
  print(c("file_path_df", file_path_df))
  #we write the general community/node overview to file
  write.csv(obj, file = paste(file_path_df, name_test, "/summary_statitics.csv", sep=""), row.names = FALSE, quote = FALSE)
  cat(paste("data frame generated: ",file_path_df, name_test, "/summary_statitics.csv", sep=""))
  #we write the general community/node overview to file
  write.table(obj_G, file = paste(file_path_df, name_test, "/communities_members-PD.csv", sep=""), sep = "\t", row.names = FALSE, quote = FALSE, na="", col.names = FALSE)
  cat(paste("\ndata frame generated:",file_path_df, name_test, "/communities_members-PD.csv", sep=""))
  #we write the general community/node overview to file
  #  write.table(obj_m, file = paste(file_path_df, datasets[i], "_g_spectral-micorarray_community.csv", sep=""), sep = "\t", row.names = FALSE, quote = FALSE, na="", col.names = FALSE)
  #  cat(paste("data frame generated:",file_path_df, datasets[i], "_g_spectral-micorarray_community.csv", sep=""))
}

#we generate a function to save community information 
community_outfile <- function(graph_biggest_component, name_test){
  #we generate a dataframe to fill with EntrezID and community information 
  obj <- data.frame(matrix(nrow=length(vertex_attr(graph_biggest_component)$name), ncol = 2))
  #we add names
  names(obj)<-c("EntrezID", "Community")
  #we add the names to the first colum and community numbers to the second
  obj$EntrezID <- vertex_attr(graph_biggest_component)$name
  obj$Community <- vertex_attr(graph_biggest_component)$community
  #we change the names of the columns
  names(obj) <- c("#clustering", "")
  #we generate the path where the outfile will be stores
  file_path_out <- paste(homedir, "/", folder, sep="")
  #we write to the output file
  write.table(obj, file = paste(file_path_out, name_test, "/communities.csv", sep=""), row.names = FALSE, quote = FALSE, sep="\t")
  cat(paste("data frame generated: ",file_path_out, name_test, "/communities.csv", sep=""))
}


#############################################################
#we add the optparse library to read command line arguments
library("optparse")

option_list = list(
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="output folder name [default= NULL]", metavar="character"),
  make_option(c("-e", "--home"), type="character", default=NULL, 
              help="home directroy name [default= NULL]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Output folder was not specified - script was stopped - define folder with -o/--out to generate graphs and data tables.n", call.=FALSE)
}

#we set the home directory - this is currently hard coded
homedir <- paste(opt$home, "/002_network_out", sep="")
print(c("homedir", homedir))
#need to hardcode on eddie
#homedir <- "/exports/eddie/scratch/s1368042/"
#we set the name of the folder of interest - based on the supplied command line input
folder = paste(opt$out, "/", sep="")
print(c(folder, "folder"))

#we reopen/define the log file
log_file <- paste(homedir, "/", folder, "log_file.txt", sep = "")
#we print to the log file - 
write("We are running the 5th file - generating networks/plots and data-overview output", file=log_file, append = TRUE)

#we generate a pdf to save the plots

file_name_plots = paste(homedir, "/", folder, "plots.pdf", sep="")
write(paste("The name of the output file for the plots is: ", file_name_plots, sep =""), file=log_file, append = TRUE)
#print(c("file_name_plots", file_name_plots))
###print to pdf
pdf(file_name_plots)

#we generate a list for the PD hits in the corresponding network
PD <- c()
#we import the protein table of interest (nodes for the network)
#it contains general information about gene-PD association
#BUT community information is specific for spectral algorithm!
prot_table = read.csv(paste(homedir, "/", folder, "node_table_igraph.csv" , sep=""), sep="\t", header=T)
#we import the edge table
edge_table = read.csv(paste(homedir, "/", folder, "PPIs_interest.csv", sep=""), sep="\t", header=F)
print("node and edge table are read")
write("node and edge table are read", file=log_file, append = TRUE)
#we generate the graph
g = graph_from_data_frame(edge_table, directed = FALSE, vertices = prot_table)

#we set vertex attributes, e.g. colour to grey
g <- set_vertex_attr(g, "color", V(g), "darkgrey")
g <- set_vertex_attr(g, "label.color", V(g), "darkgrey")
g <- set_edge_attr(g, "edge.color", E(g), "snow3")
g <- set_vertex_attr(g, "label.degree", V(g), 0)
g <- set_vertex_attr(g, "label.dist", V(g), 0.5)
#we iter over the nodes and change the colour of PD-trait associated ones
#additionally we fill the list with PD-trait genes
for(a in 1:length(vertex_attr(g)$color)){
  if(vertex_attr(g)$PD[a] == 1){
    #we add the gene name to the corresponding list
    PD <- c(PD, vertex_attr(g)$name[a])
    #we change the node colour
    vertex_attr(g)$color[a] <- "red"
    vertex_attr(g)$label.color[a] <- "red"
  }
}
#we plot the first results
plot(g, vertex.label = V(g)$gene_name, edge.color = E(g)$edge.color, vertex.size=5, vertex.label.cex = 0.7, vertex.label.color = V(g)$label.color, main=paste("full graph", sep = ""), vertex.shape= V(g)$vertex.shapes, vertex.color = V(g)$color, sub = paste("vertices: ", gorder(g), " edges: ", nrow(as_edgelist(g)), " PD-hits: ", sum(vertex_attr(g)$PD), sep = ""))
#we print to the log file
write("the first plot was printed", file=log_file, append=TRUE)

#vertex.label.degree = V(g)$label.degree, 
#we print PD hits in network
write(paste("The number of PD hits in the full network is: ", sum(vertex_attr(g)$PD), sep="") , file=log_file, append = TRUE)

#we decompose the network to define the biggest connected component
#CLUSTERING
#finding the maximal connected components (weakly or strongly)
g_components <- components(g)
#mode can be "weak" or "strong" --> current cluster - no differences 
#we decompose the graph g and select the biggest component/cluster
g_biggest_component <- decompose.graph(g)[[which.max(g_components$csize)]]

#we print the graph summary to the log file
write(paste("The number of edges in the biggest component is: ", ecount(g_biggest_component) , sep = " "), file=log_file, append=TRUE)
write(paste("The number of nodes/vertices in the biggest component is: ", vcount(g_biggest_component) , sep = " "), file=log_file, append=TRUE)

algorithm_name <- "spectral"
write(paste("The current algorithm is: ", algorithm_name , sep = " "), file=log_file, append=TRUE)
#we generate a subfolder with the algorithm name 
dir.create(paste(homedir, folder, algorithm_name, sep = "/"))
#we start with the spectral algorithm 

#set the general layout of the graph
lg <- layout_(g_biggest_component, in_circle()) * 10
#add kamada kawai
#we call the functions abofe - setting layout, generating output tables and the clustered graph
lcg <- layout_with_kk(g_biggest_component) * 25
#we call the function to set the final network layout and generate the plots
write("we are generating the network of the biggest component", file=log_file, append=TRUE)
network_layout("spectral", g_biggest_component, algorithm_name)
write("we generate dataframes for the communities", file = log_file, append = TRUE)
dataframe_communities(g_biggest_component, algorithm_name)
write("we generate a dataframe for the proteins", file = log_file, append = TRUE)
dataframe_proteins(g_biggest_component, algorithm_name)
#we generate the communities outfile
write("we re-generated the communities outfile", file = log_file, append = TRUE)
community_outfile(g_biggest_component, algorithm_name)
write(paste("all functions are finished - the data and pdf are now stored in the outfolder: ", homedir, "/", folder, sep=""), file = log_file, append = TRUE)

#we generate a dataframe to store all the nodes in the biggest connected components and the communities they belong to when applying the different 
#clustering algorithms

#we generate a matrix with one row per gene in the biggest connected component
community_dataframe <- data.frame(matrix(nrow=vcount(g_biggest_component), ncol = 1))
#we add the column name 
names(community_dataframe) <- c("GeneID")
#we add all nodes to the first column "GeneID
community_dataframe$GeneID <- get.vertex.attribute(g_biggest_component)$name
#we add the community information to the dataframe before over-writing it with community information from the next algorithm
#we add all nodes to the first column "GeneID
community_dataframe$spectral <- get.vertex.attribute(g_biggest_component)$community


##################################
#we set the algorithm name
algorithm_name <- "fast_greedy"
write(paste("The current algorithm is: ", algorithm_name , sep = " "), file=log_file, append=TRUE)
#we generate a subfolder with the algorithm name 
dir.create(paste(homedir, folder, algorithm_name, sep = "/"))

#find the communities
#Many networks consist of modules which are densely connected themselves but sparsely connected to other modules.
gbc_fast_greedy <- cluster_fast_greedy(simplify(g_biggest_component), merges = TRUE, modularity = TRUE, membership = TRUE, weights = NULL)
#we add the community based on the Newman_Girvan edge betweenness clustering algorithm to as a vertex attribute
#therefore we generate a new attribute and add the respective values
#we set the community vertex attribute 
g_biggest_component <- set_vertex_attr(g_biggest_component, "fast_greedy", V(g_biggest_component), NA)
#we transform the data into a dataframe
for(ID in community_dataframe$GeneID){
  #we print the GeneID 
  #we iter over the different communities and check the presence of the protein 
  for(e in 1:length(gbc_fast_greedy)){
    # print(c("e", e))
    if(ID %in% gbc_fast_greedy[[e]]){
      #community_dataframe$community[community_dataframe$GeneID == ID] <- e
      vertex_attr(g_biggest_component)$community[vertex_attr(g_biggest_component)$name == ID] <- e
      #print(c("e", e))
    }
  }
}

#we add the community information to the dataframe before over-writing it with community information from the next algorithm
#we add all nodes to the first column "GeneID
community_dataframe$fast_greedy <- get.vertex.attribute(g_biggest_component)$community

#we print the graph summary to the log file
write(paste("The number of edges in the biggest component is: ", ecount(g_biggest_component) , sep = " "), file=log_file, append=TRUE)
write(paste("The number of nodes/vertices in the biggest component is: ", vcount(g_biggest_component) , sep = " "), file=log_file, append=TRUE)

#set the general layout of the graph
lg <- layout_(g_biggest_component, in_circle()) * 10
#add kamada kawai
#we call the functions above - setting layout, generating output tables and the clustered graph
lcg <- layout_with_kk(g_biggest_component) * 25
#we call the function to set the final network layout and generate the plots
write("we are generating the network of the biggest component", file=log_file, append=TRUE)
#we run the network layout function with the Newman_Girvan results
network_layout(gbc_fast_greedy, g_biggest_component, algorithm_name)
write("we generate dataframes for the communities", file = log_file, append = TRUE)
dataframe_communities(g_biggest_component, algorithm_name)
write("we generate a dataframe for the proteins", file = log_file, append = TRUE)
dataframe_proteins(g_biggest_component, algorithm_name)
#we generate the communities outfile
write("we generated the communities outfile", file = log_file, append = TRUE)
community_outfile(g_biggest_component, algorithm_name)
write(paste(algorithm_name, " all functions are finished - the data and pdf are now stored in the outfolder: ", homedir, "/", folder, sep=""), file = log_file, append = TRUE)

##################################
#we set the algorithm name
algorithm_name <- "spinglas"
write(paste("The current algorithm is: ", algorithm_name , "default step size = 4", sep = " "), file=log_file, append=TRUE)
#we generate a subfolder with the algorithm name 
dir.create(paste(homedir, folder, algorithm_name, sep = "/"))

#find the communities
#Many networks consist of modules which are densely connected themselves but sparsely connected to other modules.
gbc_spinglas <-  cluster_spinglass(simplify(g_biggest_component), weights = NULL, vertex = NULL, spins = 25, parupdate = FALSE, start.temp = 1, stop.temp = 0.01, cool.fact = 0.99, update.rule = c("config", "random", "simple"), gamma = 1, implementation = c("orig", "neg"), gamma.minus = 1)
#we add the community based on the Newman_Girvan edge betweenness clustering algorithm to as a vertex attribute
#therefore we generate a new attribute and add the respective values
#we set the community vertex attribute 
g_biggest_component <- set_vertex_attr(g_biggest_component, "spinglas", V(g_biggest_component), NA)
#we transform the data into a dataframe
for(ID in community_dataframe$GeneID){
  #we print the GeneID 
  #we iter over the different communities and check the presence of the protein 
  for(e in 1:length(gbc_spinglas)){
    # print(c("e", e))
    if(ID %in% gbc_spinglas[[e]]){
      #community_dataframe$community[community_dataframe$GeneID == ID] <- e
      vertex_attr(g_biggest_component)$community[vertex_attr(g_biggest_component)$name == ID] <- e
      #print(c("e", e))
    }
  }
}

#we add the community information to the dataframe before over-writing it with community information from the next algorithm
#we add all nodes to the first column "GeneID
community_dataframe$spinglas <- get.vertex.attribute(g_biggest_component)$community

#we print the graph summary to the log file
write(paste("The number of edges in the biggest component is: ", ecount(g_biggest_component) , sep = " "), file=log_file, append=TRUE)
write(paste("The number of nodes/vertices in the biggest component is: ", vcount(g_biggest_component) , sep = " "), file=log_file, append=TRUE)

#set the general layout of the graph
lg <- layout_(g_biggest_component, in_circle()) * 10
#add kamada kawai
#we call the functions above - setting layout, generating output tables and the clustered graph
lcg <- layout_with_kk(g_biggest_component) * 25
#we call the function to set the final network layout and generate the plots
write("we are generating the network of the biggest component", file=log_file, append=TRUE)
#we run the network layout function with the Newman_Girvan results
network_layout(gbc_spinglas, g_biggest_component, algorithm_name)
write("we generate dataframes for the communities", file = log_file, append = TRUE)
dataframe_communities(g_biggest_component, algorithm_name)
write("we generate a dataframe for the proteins", file = log_file, append = TRUE)
dataframe_proteins(g_biggest_component, algorithm_name)
#we generate the communities outfile
write("we generated the communities outfile", file = log_file, append = TRUE)
community_outfile(g_biggest_component, algorithm_name)
write(paste(algorithm_name, " all functions are finished - the data and pdf are now stored in the outfolder: ", homedir, "/", folder, sep=""), file = log_file, append = TRUE)


##################################
#we set the algorithm name
algorithm_name <- "louvain"
write(paste("The current algorithm is: ", algorithm_name, sep = " "), file=log_file, append=TRUE)
#we generate a subfolder with the algorithm name 
dir.create(paste(homedir, folder, algorithm_name, sep = "/"))

#find the communities
#Many networks consist of modules which are densely connected themselves but sparsely connected to other modules.
#louvain algorithm - multi-level modularity optimization algorithm
gbc_louvain <- cluster_louvain(simplify(g_biggest_component), weights = NULL)
#we add the community based on the Newman_Girvan edge betweenness clustering algorithm to as a vertex attribute
#therefore we generate a new attribute and add the respective values
#we set the community vertex attribute 
g_biggest_component <- set_vertex_attr(g_biggest_component, "spinglas", V(g_biggest_component), NA)
#we transform the data into a dataframe
for(ID in community_dataframe$GeneID){
  #we print the GeneID 
  #we iter over the different communities and check the presence of the protein 
  for(e in 1:length(gbc_louvain)){
    # print(c("e", e))
    if(ID %in% gbc_louvain[[e]]){
      #community_dataframe$community[community_dataframe$GeneID == ID] <- e
      vertex_attr(g_biggest_component)$community[vertex_attr(g_biggest_component)$name == ID] <- e
      #print(c("e", e))
    }
  }
}

#we add the community information to the dataframe before over-writing it with community information from the next algorithm
#we add all nodes to the first column "GeneID
community_dataframe$louvain <- get.vertex.attribute(g_biggest_component)$community

#we print the graph summary to the log file
write(paste("The number of edges in the biggest component is: ", ecount(g_biggest_component) , sep = " "), file=log_file, append=TRUE)
write(paste("The number of nodes/vertices in the biggest component is: ", vcount(g_biggest_component) , sep = " "), file=log_file, append=TRUE)

#set the general layout of the graph
lg <- layout_(g_biggest_component, in_circle()) * 10
#add kamada kawai
#we call the functions above - setting layout, generating output tables and the clustered graph
lcg <- layout_with_kk(g_biggest_component) * 25
#we call the function to set the final network layout and generate the plots
write("we are generating the network of the biggest component", file=log_file, append=TRUE)
#we run the network layout function with the Newman_Girvan results
network_layout(gbc_louvain, g_biggest_component, algorithm_name)
write("we generate dataframes for the communities", file = log_file, append = TRUE)
dataframe_communities(g_biggest_component, algorithm_name)
write("we generate a dataframe for the proteins", file = log_file, append = TRUE)
dataframe_proteins(g_biggest_component, algorithm_name)
#we generate the communities outfile
write("we generated the communities outfile", file = log_file, append = TRUE)
community_outfile(g_biggest_component, algorithm_name)
write(paste(algorithm_name, " all functions are finished - the data and pdf are now stored in the outfolder: ", homedir, "/", folder, sep=""), file = log_file, append = TRUE)



##################################
#we set the algorithm name
algorithm_name <- "infomap"
write(paste("The current algorithm is: ", algorithm_name, sep = " "), file=log_file, append=TRUE)
#we generate a subfolder with the algorithm name 
dir.create(paste(homedir, folder, algorithm_name, sep = "/"))

#find the communities
#Many networks consist of modules which are densely connected themselves but sparsely connected to other modules.
#infomap algorithm - community structure minimizes the expected description length of a random walker trajectory
gbc_infomap <- cluster_infomap(simplify(g_biggest_component), e.weights = NULL, v.weights = NULL, nb.trials = 10, modularity = TRUE)
#we add the community based on the Newman_Girvan edge betweenness clustering algorithm to as a vertex attribute
#therefore we generate a new attribute and add the respective values
#we set the community vertex attribute 
g_biggest_component <- set_vertex_attr(g_biggest_component, "infomap", V(g_biggest_component), NA)
#we transform the data into a dataframe
for(ID in community_dataframe$GeneID){
  #we print the GeneID 
  #we iter over the different communities and check the presence of the protein 
  for(e in 1:length(gbc_infomap)){
    # print(c("e", e))
    if(ID %in% gbc_infomap[[e]]){
      #community_dataframe$community[community_dataframe$GeneID == ID] <- e
      vertex_attr(g_biggest_component)$community[vertex_attr(g_biggest_component)$name == ID] <- e
      #print(c("e", e))
    }
  }
}

#we add the community information to the dataframe before over-writing it with community information from the next algorithm
#we add all nodes to the first column "GeneID
community_dataframe$infomap <- get.vertex.attribute(g_biggest_component)$community

#we print the graph summary to the log file
write(paste("The number of edges in the biggest component is: ", ecount(g_biggest_component) , sep = " "), file=log_file, append=TRUE)
write(paste("The number of nodes/vertices in the biggest component is: ", vcount(g_biggest_component) , sep = " "), file=log_file, append=TRUE)

#set the general layout of the graph
lg <- layout_(g_biggest_component, in_circle()) * 10
#add kamada kawai
#we call the functions above - setting layout, generating output tables and the clustered graph
lcg <- layout_with_kk(g_biggest_component) * 25
#we call the function to set the final network layout and generate the plots
write("we are generating the network of the biggest component", file=log_file, append=TRUE)
#we run the network layout function with the Newman_Girvan results
network_layout(gbc_infomap, g_biggest_component, algorithm_name)
write("we generate dataframes for the communities", file = log_file, append = TRUE)
dataframe_communities(g_biggest_component, algorithm_name)
write("we generate a dataframe for the proteins", file = log_file, append = TRUE)
dataframe_proteins(g_biggest_component, algorithm_name)
#we generate the communities outfile
write("we generated the communities outfile", file = log_file, append = TRUE)
community_outfile(g_biggest_component, algorithm_name)
write(paste(algorithm_name, " all functions are finished - the data and pdf are now stored in the outfolder: ", homedir, "/", folder, sep=""), file = log_file, append = TRUE)


##################################
##################################
#we close and print the plot 
dev.off()

#we also print the dataframe with node-ID and community membership for all of the algorithms
write.table(community_dataframe, file =paste(homedir, "/", folder, "community_membership_overview.csv", sep = ""), sep="\t", row.names = FALSE)

