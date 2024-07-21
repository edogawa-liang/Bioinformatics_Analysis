## ----Part1
<<<<<<< Updated upstream
genes <- read.csv("Analysis/Data/genes.csv")
Response <- read.csv("Analysis/Data/Response.csv")
=======
library(tidyverse)
library(Rtsne)
genes <- read.csv("Data/genes.csv")
metabolites <- read.csv("Data/metabolites.csv")
Response <- read.csv("Data/Response.csv")
>>>>>>> Stashed changes
prim_disease <- Response$primary_disease
table(prim_disease)
nrow(unique(Response))
nrow(unique(metabolites))
length(unique(genes$X))
ID <- unique(Response$ID)
genes_filter <- genes %>% filter(X %in% ID)
Response_filter <- unique(Response)
metabolites_filter <- unique(metabolites)
## ---- T-sne
t_genes <- t(genes_filter)
t_genes <- t_genes[2:nrow(t_genes),]
t_genes <- apply(t_genes,2 , as.numeric)
t_genes <- data.frame(t_genes)
tsne_out <- Rtsne(t_genes,dims = 3,check_duplicates = FALSE)
tsne_plot <- data.frame(x = tsne_out$Y[,1], 
                        y = tsne_out$Y[,2],
                        z = tsne_out$Y[,3])
library(plotly)
plot_ly(x=tsne_plot$x,y=tsne_plot$y,z=tsne_plot$z,type = "scatter3d", mode="markers")
# Plotting the plot using ggplot() function
ggplot2::ggplot(tsne_plot,label=Species) + geom_point(aes(x=x,y=y),color = kmedoid.cluster$clustering)

ggplot2::ggplot(tsne_plot,label=Species) + geom_point(aes(x=x,y=z),color = km$cluster)

ggplot2::ggplot(tsne_plot,label=Species) + geom_point(aes(x=z,y=y),color = km$cluster)

library(cluster)
kmedoid.cluster <- pam(t_genes, k=5)
kmedoid.cluster$clustering
kmedoid.cluster$objective

km <- kmeans(t_genes,centers = 6,nstart = 10)
km$cluster
## ----Response dealing
bone <- ifelse(Response$bone=="NonMetastasis",0,1)
brain <- ifelse(Response$brain=="NonMetastasis",0,1)
kidney <- ifelse(Response$kidney=="NonMetastasis",0,1)
liver <- ifelse(Response$liver=="NonMetastasis",0,1)
lung <- ifelse(Response$lung=="NonMetastasis",0,1)
response0_1 <- data.frame(bone,brain,kidney,liver,lung)
response_label <- c()
for (i in 1:nrow(response0_1)) {
  rr <- response0_1[i,]
  label <- rr$bone*(2^0) + rr$brain*(2^1) + rr$kidney*(2^2) + rr$liver*(2^3) + rr$lung*(2^4)
  response_label <- c(label,response_label)
}

table(response_label)