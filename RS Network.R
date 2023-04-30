library(data.table)
library(igraph)
library(data.table)
library(GGally)
library(ggplot2)
library(network)
library(sna)

###    D' Analysis    ###

#Reading files in each chromosome

D <- c("chr1D.csv","chr2D.csv","chr3D.csv",
       "chr4D.csv","chr5D.csv","chr6D.csv",
       "chr7D.csv","chr8D.csv","chr9D.csv",
       "chr10D.csv","chr11D.csv","chr12D.csv",
       "chr13D.csv","chr14D.csv","chr15D.csv",
       "chr16D.csv","chr17D.csv","chr18D.csv",
       "chr19D.csv","chr20D.csv","chr21D.csv",
       "chr22D.csv","chrXD.csv")

resultD <- list()

#Purring D' data

for (n in 1:length(D)) {
  chrD <- read.csv(D[n],row.names = 1,header = TRUE)
  len <- length(chrD)
  for (i in 1:len) {
    for (j in 1:len) {
      if (chrD[i,j] < 0.5){
        chrD[i,j] = 0}
    }
  }
  
  resultD[[n]] <- chrD
  
}



###    R2 Analysis    ###

#Reading files in each chromosome

R2 <- c("chr1R.csv","chr2R.csv","chr3R.csv",
        "chr4R.csv","chr5R.csv","chr6R.csv",
        "chr7R.csv","chr8R.csv","chr9R.csv",
        "chr10R.csv","chr11R.csv","chr12R.csv",
        "chr13R.csv","chr14R.csv","chr15R.csv",
        "chr16R.csv","chr17R.csv","chr18R.csv",
        "chr19R.csv","chr20R.csv","chr21R.csv",
        "chr22R.csv","chrXR.csv")

resultR2 <- list()

#Purring R2 data

for (n in 1:length(R2)) {
  chrR2 <- read.csv(R2[n],row.names = 1,header = TRUE)
  len <- length(chrR2)
  for (i in 1:len) {
    for (j in 1:len) {
      if (chrR2[i,j] < 0.16){
        chrR2[i,j] = 0}
    }
  }
  
  resultR2[[n]] <- chrR2
  
}



###    Network Analysis    ###

#Data prepration

result <- list()

#Reading files in each chromosome
for (n in 1:length(resultD)) {
  chrD <- resultD[[n]]
  chrR2 <- resultR2[[n]]
  len <- length(chrD)
  for (i in 1:len) {
    for (j in 1:len) {
      if (chrD[i,j] != 0 & chrR2[i,j] == 0){
        chrD[i,j] = 0}
    }
  }
  
  result[[n]] <- chrD
  
}



#Removing same cells

resultRem <- list()

for (n in 1:length(resultD)) {
  chr <- result[[n]]
  len <- length(chr)
  s <- 1
  x <- 2
  while(x <= len){
    for (i in x:len) {
      chr[i,s] <- 0
    }
    s <- s+1
    x <- x+1
  }
  resultRem[[n]] <- chr

}  


#Removing Zero rows and columns

resultZero <- list()

for (n in 1:length(resultRem)) {
  chr <- resultRem[[n]]
  len <- length(chr)
  for (i in 1:len) {
    for (j in 1:len) {
      if (chr[i,j] == 1 & i == j){ 
        chr[i,j] <- 0
      }
    }
  }
  p <- c()
  s <- 1
  x <- 1
  t <- 0
  while(x <= len){
    for (j in 1:len) {
      if (chr[x,j] == 0){
        t <- t+1
      }
    }
    if (t==len){
      p[s] <- x
      s <- s+1
    }
    x <- x+1
    t <- 0
  } 
  
  b <- c()
  s <- 1
  x <- 1
  t <- 0
  while(x <= len){
    for (j in 1:len) {
      if (chr[j,x] == 0){
        t <- t+1
      }
    }
    if (t==len){
      b[s] <- x
      s <- s+1
    }
    x <- x+1
    t <- 0
  }  
  
  chr <- chr[-p,-b]
  resultZero[[n]] <- chr
}


#Network prepration

resultNet <- list()

for (n in 1:length(resultZero)) {
  chr <- resultZero[[n]]
  len <- length(chr)
  k <- 1
  RSnetwork <- matrix(data=NA,nrow = 100000,ncol = 3,byrow = FALSE)
  for (i in 1:nrow(chr)) {
    for (j in 1:ncol(chr)) {
      if (chr[i,j] != 0){
        RSnetwork[k,1] <- rownames(chr[i,])
        RSnetwork[k,2] <- colnames(chr[j])
        RSnetwork[k,3] <- chr[i,j]
        k <- k+1
      }
    
  }
  
  RSnetworkf <- na.omit(RSnetwork)
  resultNet[[n]] <- RSnetworkf
  }
}

write.table(resultNet, file = "RSNet.txt",row.names = FALSE,col.names =FALSE, sep = "\t")



#Network construction

My_data <- resultNet
My_metrix <- as.matrix(My_data)
My_network <- graph.adjacency(My_metrix,mode = "undirected" ,weighted=TRUE, diag = FALSE)

#plot(My_network,layout = layout.circle)
ggnet2(My_network)
plot(My_network,layout = layout.fruchterman.reingold, vertex.color = "cyan")
plot(My_network,layout = layout.graphopt, vertex.color = "cyan")
plot(My_network , layout = layout.kamada.kawai ,vertex.color = "cyan")
edge.betweenness(My_network)


plot(My_network , vertex.color = "cyan")

png(filename = "rsNetwork.png")
dev.off()
max_cliques(My_network, min = 4, max = NULL, subset = NULL, file = NULL)


hub_score(My_network, scale = TRUE, weights = NULL,
          options = arpack_defaults)
max(degree(My_network, v = V(My_network), mode = c("all", "out", "in", "total"),
           loops = TRUE, normalized = FALSE))
sum(degree_distribution(My_network, cumulative = FALSE))
E(My_network)
V(My_network)
ecount(My_network)
vcount(My_network)

