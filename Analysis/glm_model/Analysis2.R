## preprocess
t_metabolites_filter <- metabolites_filter %>% 
  arrange(ID) %>%
  `rownames<-`(ID) %>%
  t %>%
  as.data.frame %>% 
  dplyr::filter(!row_number() %in% c(1)) %>%
  apply(2 , as.numeric)

## dimesion reducition plot
metabolites_tsne2 <- Rtsne::Rtsne(t_metabolites_filter,dims = 2,check_duplicates = FALSE)
metabolites_tsne_2data <- data.frame(x = metabolites_tsne2$Y[,1], 
                                    y = metabolites_tsne2$Y[,2])
ggplot2::ggplot(metabolites_tsne_2data,label=Species) + geom_point(aes(x=x,y=y))

metabolites_tsne3 <- Rtsne::Rtsne(t_metabolites_filter,dims = 3,check_duplicates = FALSE)
metabolites_tsne_3data <- data.frame(x = metabolites_tsne3$Y[,1], 
                                     y = metabolites_tsne3$Y[,2],
                                     z = metabolites_tsne3$Y[,3])
plotly::plot_ly(x=metabolites_tsne_3data$x,y=metabolites_tsne_3data$y,z=metabolites_tsne_3data$z,type = "scatter3d", mode="markers")

## clustering
## do clustering and see how it performs in 3 dim plot
library(cluster)
kmedoid.cluster5 <- pam(t_metabolites_filter, k=5)
kmedoid.cluster6 <- pam(t_metabolites_filter, k=6)
kmedoid.cluster7 <- pam(t_metabolites_filter, k=7)
kmedoid.cluster8 <- pam(t_metabolites_filter, k=8)
plotly::plot_ly(x=metabolites_tsne_3data$x,y=metabolites_tsne_3data$y,z=metabolites_tsne_3data$z,type = "scatter3d", mode="markers",color = kmedoid.cluster5$clustering)
plotly::plot_ly(x=metabolites_tsne_3data$x,y=metabolites_tsne_3data$y,z=metabolites_tsne_3data$z,type = "scatter3d", mode="markers",color = kmedoid.cluster6$clustering)
plotly::plot_ly(x=metabolites_tsne_3data$x,y=metabolites_tsne_3data$y,z=metabolites_tsne_3data$z,type = "scatter3d", mode="markers",color = kmedoid.cluster7$clustering)
plotly::plot_ly(x=metabolites_tsne_3data$x,y=metabolites_tsne_3data$y,z=metabolites_tsne_3data$z,type = "scatter3d", mode="markers",color = kmedoid.cluster8$clustering)
### using kmeans
km6 <- kmeans(t_metabolites_filter,centers = 6,nstart = 15)
plotly::plot_ly(x=metabolites_tsne_3data$x,y=metabolites_tsne_3data$y,z=metabolites_tsne_3data$z,type = "scatter3d", mode="markers",color = km6$cluster)
km5 <- kmeans(t_metabolites_filter,centers = 5,nstart = 15)
plotly::plot_ly(x=metabolites_tsne_3data$x,y=metabolites_tsne_3data$y,z=metabolites_tsne_3data$z,type = "scatter3d", mode="markers",color = km5$cluster)
km4 <- kmeans(t_metabolites_filter,centers = 4,nstart = 15)
plotly::plot_ly(x=metabolites_tsne_3data$x,y=metabolites_tsne_3data$y,z=metabolites_tsne_3data$z,type = "scatter3d", mode="markers",color = km4$cluster)

### Define metobolites cluster
## By clustering algorithm, we can obtain the information that which metabolites have similar information w.r.t all cells.
## We'll define a function that groups metabolites.
### According to the groups of metabolites, we have to reduce it to one dimension.
### In our case, we will take medium of all metabolites of a cell to represent the group
## Reduce by mean
reduce_metabolites.df <- t_metabolites_filter %>% 
  data.frame %>%
  mutate(cluster  = km5$cluster, metabolites = colnames(metabolites_filter)[2:ncol(metabolites_filter)]) %>%
  group_by(cluster) %>%
  summarise(across(where(is.numeric), mean))

### analysis data
## 1. primary disease convert
Response_filter <- Response_filter %>%
  arrange(ID)
prim_dis <- ifelse(Response_filter$primary_disease %in% names(which(table(Response_filter$primary_disease)<5)),"Other disease",Response_filter$primary_disease)
## 1.1 convert new category
bone <- ifelse(Response_filter$bone=="NonMetastasis",0,1)
brain <- ifelse(Response_filter$brain=="NonMetastasis",0,1)
kidney <- ifelse(Response_filter$kidney=="NonMetastasis",0,1)
liver <- ifelse(Response_filter$liver=="NonMetastasis",0,1)
lung <- ifelse(Response_filter$lung=="NonMetastasis",0,1)
response0_1 <- data.frame(bone,brain,kidney,liver,lung)
response_label <- c()
for (i in 1:nrow(response0_1)) {
  rr <- response0_1[i,]
  label <- rr$bone*(2^0) + rr$brain*(2^1) + rr$kidney*(2^2) + rr$kidney*(2^3) + rr$lung*(2^4)
  response_label <- c(label,response_label)
}

table(response_label)
## 2. data.frame
df.analysis <- t(reduce_metabolites.df)[2:ncol(reduce_metabolites.df),] %>%
  data.frame %>%
  mutate(bone = as.factor(Response_filter$bone),brain = as.factor(Response_filter$brain),kidney = as.factor(Response_filter$kidney),liver = as.factor(Response_filter$liver),lung = as.factor(Response_filter$lung),primary_disease = as.factor(prim_dis),response_label = as.factor(response_label))

## Build logistic regression
glimpse(df.analysis)
table(df.analysis$primary_disease)
attach(df.analysis)
logit.bone <- glm(bone ~ X1+X2+X3+X4+X5+primary_disease,data = df.analysis,family = "binomial")
summary(logit.bone)

logit.brain <- glm(brain ~ X1+X2+X3+X4+X5+primary_disease,data = df.analysis,family = "binomial")
summary(logit.brain)

logit.kidney <- glm(kidney ~ X1+X2+X3+X4+X5+primary_disease,data = df.analysis,family = "binomial")
summary(logit.kidney)

logit.liver <- glm(liver ~ X1+X2+X3+X4+X5+primary_disease,data = df.analysis,family = "binomial")
summary(logit.liver)

logit.lung <- glm(lung ~ X1+X2+X3+X4+X5+primary_disease,data = df.analysis,family = "binomial")
summary(logit.lung)


### Analysis 2
library(nnet)
attach(df.analysis)
multilogit <- multinom(response_label ~ X1+X2+X3+X4+X5+primary_disease,data = df.analysis)
summary(multilogit)






## zero inflated model
## 可以先預測轉移到幾處，在處理轉移到哪裡
## 資料型態會是轉移到幾處，所以先把資料變成count資料
bone <- ifelse(Response_filter$bone=="NonMetastasis",0,1)
brain <- ifelse(Response_filter$brain=="NonMetastasis",0,1)
kidney <- ifelse(Response_filter$kidney=="NonMetastasis",0,1)
liver <- ifelse(Response_filter$liver=="NonMetastasis",0,1)
lung <- ifelse(Response_filter$lung=="NonMetastasis",0,1)
response0_1 <- data.frame(bone,brain,kidney,liver,lung)
response_label <- c()
for (i in 1:nrow(response0_1)) {
  rr <- response0_1[i,]
  label <- rr$bone + rr$brain + rr$kidney + rr$kidney + rr$lung
  response_label <- c(label,response_label)
}

table(response_label)
df.analysis <- t(reduce_metabolites.df)[2:ncol(reduce_metabolites.df),] %>%
  data.frame %>%
  mutate(bone = as.factor(Response_filter$bone),brain = as.factor(Response_filter$brain),kidney = as.factor(Response_filter$kidney),liver = as.factor(Response_filter$liver),lung = as.factor(Response_filter$lung),primary_disease = as.factor(prim_dis),response_label = as.integer(response_label))
library(pscl)
glimpse(df.analysis)
attach(df.analysis)
zimodel <- zeroinfl(response_label~X1+X2+X3+X4+X5|X1+X2+X3+X4+X5,data = df.analysis)
summary(zimodel)
hurdlemodel <- hurdle(response_label~X1+X2+X3+X4+X5|X1+X2+X3+X4+X5,data = df.analysis)
summary(hurdlemodel)


