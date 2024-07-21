#setwd("Analysis")
genes <- read.csv("Data/genes.csv")
metabolites <- read.csv("Data/metabolites.csv")
Response <- read.csv("Data/Response.csv")
prim_disease <- Response$primary_disease
table(prim_disease)
nrow(unique(Response))
nrow(unique(metabolites))
length(unique(genes$X))
ID <- unique(Response$ID)
genes_filter <- genes %>% 
  filter(X %in% ID) %>% 
  arrange(X)
Response_filter <- unique(Response) %>%
  arrange(ID)
metabolites_filter <- unique(metabolites) %>%
  arrange(ID) %>%
  select(-c(ID))

thresh <- colSums(genes_filter > 0.5) >= 2
genes_cut <- genes_filter[,thresh]
genes_cut <- genes_cut[,2:ncol(genes_cut)]
gene_sd <- apply(genes_cut, 2, sd)
extreme_thresh <- abs(t(genes_cut) - colMeans(genes_cut)) > (3 * gene_sd)
extreme_thresh <- colSums(t(extreme_thresh)) > 0
genes_cut <- genes_cut[,!extreme_thresh]

genes_filter <- genes_cut
save(genes_filter,Response_filter,metabolites_filter,file = "Datas.RData")

## ----Read Data
load("Datas.RData")
##################################################
## 1. 區分出是否轉移
### 1.1 挑選出變數
## 2. 區分出轉移到幾個地方
## 3. 區分出轉移到哪個地方
##################################################
  
### 取得是否轉移的類別
# organ.vec <- c("bone","brain","kidney","liver","lung")
# for(organ in organ.vec){
#   organ.name <- as.name(organ)
#   trans.vec <- Response_filter %>% 
#     mutate(organ.name = replace(organ.name, organ.name=="NonMetastasis",0))
#   trans.vec <- ifelse(trans.vec=="NonMetastasis",0,1)
#   
#   assign(organ,)
# }

## ----Preprocess
bone <- ifelse(Response_filter$bone=="NonMetastasis",0,1)
brain <- ifelse(Response_filter$brain=="NonMetastasis",0,1)
kidney <- ifelse(Response_filter$kidney=="NonMetastasis",0,1)
liver <- ifelse(Response_filter$liver=="NonMetastasis",0,1)
lung <- ifelse(Response_filter$lung=="NonMetastasis",0,1)
response0_1 <- data.frame(bone,brain,kidney,liver,lung)
count_response_label <- c()
for (i in 1:nrow(response0_1)) {
  rr <- response0_1[i,]
  label <- rr$bone + rr$brain + rr$kidney + rr$kidney + rr$lung
  count_response_label <- c(count_response_label,label)
}
binary_response_label <- ifelse(count_response_label>=1,1,0)

## ----lasso

## 1.2 lasso for binary
library(glmnet)
set.seed(1000)
metabolites.cv_lasso_binary <- cv.glmnet(as.matrix(metabolites_filter),binary_response_label, family = "binomial", alpha = 1)

metabolites.lasso_logistic_binary <- glmnet(as.matrix(metabolites_filter),binary_response_label, family = "binomial", alpha = 1,lambda = metabolites.cv_lasso_binary$lambda.min)

metabolites.coefficients_matrix_binary <- as.matrix(coef(metabolites.lasso_logistic_binary))
selected_variables_binary <- rownames(coefficients_matrix_binary)[coefficients_matrix_binary != 0]
metabolites.selected_variables_binary


## 1.3 lasso for count
set.seed(1000)
cv_lasso_count <- cv.glmnet(as.matrix(metabolites_filter),count_response_label, family = "poisson", alpha = 1)

lasso_logistic_count <- glmnet(as.matrix(metabolites_filter),count_response_label, family = "poisson", alpha = 1,lambda = cv_lasso_count$lambda.min)

coefficients_matrix_count <- as.matrix(coef(lasso_logistic_count))
selected_variables_count <- rownames(coefficients_matrix_count)[coefficients_matrix_count != 0]
selected_variables_count

selected_variables_binary %in% selected_variables_count
selected_variables_count %in% selected_variables_binary


## ----Deal with interaction
## 1.4 create interaction dataframe
interaction.trasfer.function <- function(df){
  ## this function inputs a dataframe
  ## this function returns a dataframe with interaction term inside the dataframe
  l <- colnames(df)
  name.vec <- c()
  for(i in 1:length(l)){
    for(j in 1:length(l)){
      if(i<=j) next
      name1 <- l[i];name2 <- l[j]
      name = paste(name1,"_",name2)
      name.vec <- c(name.vec,name)
      assign(name,df[,i]*df[,j])
      df <- cbind(df,get(name))
      rm(list = c(name))
      print(length(name.vec))
    }
  }
  colnames(df) <- c(l,name.vec)
  return(df)
}
#save(interaction.trasfer.function,file = "interaction_transfer.RData")

### 2. zero inflated model
load("binary.Rdata")
load("count.Rdata")
load("response.Rdata")
load("variable.Rdata")


library(pscl)
add_symbol <- function(string){
  out_string <- paste("`",string,"`",sep = "")
  return(out_string)
}
selected_variables_binary <- sapply(variables_binary,add_symbol)
selected_variables_count <- sapply(variables_count,add_symbol)
  

X <- cbind(interaction.trasfer.function(x_binary),interaction.trasfer.function(x_count))
filter_colnames <- names(table(colnames(X)))
X <- X[,filter_colnames]

df.analysis <- cbind(X,binary_response_label,count_response_label)

print(paste(selected_variables_count[2:length(selected_variables_count)],collapse  = "+"))

zimodel <- zeroinfl(as.formula(paste("count_response_label~",paste(selected_variables_count[2:length(selected_variables_count)],collapse  = "+"),"|",paste(selected_variables_binary[2:length(selected_variables_binary)],collapse = "+"))),data = df.analysis)
summary(zimodel)
library(stargazer)
hurdlemodel <- hurdle(as.formula(paste("count_response_label~",paste(selected_variables_count[2:length(selected_variables_count)],collapse  = "+"),"|",paste(selected_variables_binary[2:length(selected_variables_binary)],collapse = "+"))),data = df.analysis)
summary(hurdlemodel)

purepossionmodel <- glm(as.formula(paste("count_response_label~",paste(selected_variables_count[2:length(selected_variables_count)],collapse  = "+"))),data = df.analysis,family = "poisson")
summary(purepossionmodel)
VIF(purepossionmodel)
stargazer(hurdlemodel,purepossionmodel,purelogitmodel)

purelogitmodel <- glm(as.formula(paste("count_response_label~",paste(selected_variables_binary[2:length(selected_variables_binary)],collapse = "+"))),data = df.analysis,family = "poisson")
summary(purelogitmodel)
VIF(purelogitmodel)

p <- ifelse(round(predict(zimodel))>5,5,round(predict(zimodel)))
mean((p-count_response_label)^2)

p <- ifelse(round(predict(hurdlemodel))>5,5,round(predict(hurdlemodel)))
mean((p-count_response_label)^2)

p <- ifelse(round(predict(purepossionmodel))>5,5,round(predict(purepossionmodel)))
mean((p-count_response_label)^2)

#### other method to pick variables
load("filter_data2.RData")
organ_candidate <- data.frame(organ_candidate) %>%
  mutate(binary_response_label,count_response_label)
selected_variables <- sapply(colnames(organ_candidate)[1:(length(colnames(organ_candidate))-2)],add_symbol)

zimodel <- zeroinfl(as.formula(paste("count_response_label~",paste(selected_variables,collapse  = "+"),"|",paste(selected_variables,collapse = "+"))),data = organ_candidate)
summary(zimodel)

hurdlemodel <- hurdle(as.formula(paste("count_response_label~",paste(selected_variables,collapse  = "+"),"|",paste(selected_variables,collapse = "+"))),data = organ_candidate)
summary(hurdlemodel)




### by 器官
## bone
bone.df <- x_bone %>%
  interaction.trasfer.function %>%
  select(all_of(variables_bone[2:length(variables_bone)])) %>%
  mutate(bone)

## brain
brain.df <- x_brain %>%
  interaction.trasfer.function %>%
  select(all_of(variables_brain[2:length(variables_brain)])) %>%
  mutate(brain)

## lung
lung.df <- x_lung %>%
  interaction.trasfer.function %>%
  select(all_of(variables_lung[2:length(variables_lung)])) %>%
  mutate(lung)

## liver
liver.df <- x_liver %>%
  interaction.trasfer.function %>%
  select(all_of(variables_liver[2:length(variables_liver)])) %>%
  mutate(liver)

## kidney
kidney.df <- x_kidney %>%
  interaction.trasfer.function %>%
  select(all_of(variables_kidney[2:length(variables_kidney)])) %>%
  mutate(kidney)

#### bone
set.seed(2000)
train.bone <- dplyr::slice_sample(bone.df,prop = 0.8,by = bone)
test.bone <- dplyr::anti_join(bone.df,train.bone)
table(train.bone$bone)
table(test.bone$bone)
logit.train.bone <- glm(bone~.,data = train.bone,family = 'binomial')
pred.test.bone <- predict(logit.train.bone,type = 'response',newdata = test.bone )
pred.bone <- ifelse(pred.test.bone>0.5,1,0)
bone_table <- table(predicted = pred.bone,actual = test.bone$bone)
bone_table
bone.pred <- predict(logit.train.bone,type = "response")
proc.bone <- roc(train.bone$bone,bone.pred)
proc.bone$auc
ggroc(proc.bone)


logit.bone <- glm(bone~.,data = bone.df,family = "binomial")
summary(logit.bone)
bone.pred <- predict(logit.bone,type = "response")
proc.bone <- roc(bone,bone.pred)
proc.bone$auc
ggroc(proc.bone)
VIF(logit.bone)


##### brain
set.seed(2000)
train.brain <- dplyr::slice_sample(brain.df,prop = 0.8,by = brain)
test.brain <- dplyr::anti_join(brain.df,train.brain)
table(train.brain$brain)
table(test.brain$brain)
logit.train.brain <- glm(brain~.,data = train.brain,family = 'binomial')
pred.test.brain <- predict(logit.train.brain,type = 'response',newdata = test.brain )
pred.brain <- ifelse(pred.test.brain>0.5,1,0)
brain_table <- table(predicted = pred.brain,actual = test.brain$brain)
brain_table
brain.pred <- predict(logit.train.brain,type = "response")
proc.brain <- roc(train.brain$brain,brain.pred)
proc.brain$auc
ggroc(proc.brain)


logit.brain <- glm(brain~.,data = brain.df,family = "binomial")
summary(logit.brain)
brain.pred <- predict(logit.brain,type = "response")
proc.brain <- roc(brain,brain.pred)
proc.brain$auc
ggroc(proc.brain)

#### lung
set.seed(2000)
train.lung <- dplyr::slice_sample(lung.df,prop = 0.8,by = lung)
test.lung <- dplyr::anti_join(lung.df,train.lung)
table(train.lung$lung)
table(test.lung$lung)
logit.train.lung <- glm(lung~.,data = train.lung,family = 'binomial')
pred.test.lung <- predict(logit.train.lung,type = 'response',newdata = test.lung )
pred.lung <- ifelse(pred.test.lung>0.5,1,0)
lung_table <- table(predicted = pred.lung,actual = test.lung$lung)
lung_table
lung.pred <- predict(logit.train.lung,type = "response")
proc.lung <- roc(train.lung$lung,lung.pred)
proc.lung$auc
ggroc(proc.lung)


logit.lung <- glm(lung~.,data = lung.df,family = "binomial")
summary(logit.lung)
lung.pred <- predict(logit.lung,type = "response")
proc.lung <- roc(lung,lung.pred)
proc.lung$auc
ggroc(proc.lung)

#### liver
logit.liver <- glm(liver~.,data = liver.df,family = "binomial",maxit=1000)
summary(logit.liver)
liver.pred <- predict(logit.liver,type = "response")
proc.liver <- roc(liver,liver.pred)
proc.liver$auc
ggroc(proc.liver)


#### kidney
set.seed(2000)
train.kidney <- dplyr::slice_sample(kidney.df,prop = 0.8,by = kidney)
test.kidney <- dplyr::anti_join(kidney.df,train.kidney)
table(train.kidney$kidney)
table(test.kidney$kidney)
logit.train.kidney <- glm(kidney~.,data = train.kidney,family = 'binomial',maxit = 3000)
pred.test.kidney <- predict(logit.train.kidney,type = 'response',newdata = test.kidney )
pred.kidney <- ifelse(pred.test.kidney>0.5,1,0)
kidney_table <- table(predicted = pred.kidney,actual = test.kidney$kidney)
aaaa <- roc(test.kidney$kidney,pred.test.kidney)
aaaa$auc

kidney_table
kidney.pred <- predict(logit.train.kidney,type = "response")
proc.kidney <- roc(train.kidney$kidney,kidney.pred)
proc.kidney$auc
ggroc(proc.kidney)


logit.kidney <- glm(kidney~.,data = kidney.df,family = "binomial")
summary(logit.kidney)
kidney.pred <- predict(logit.kidney,type = "response")
proc.kidney <- roc(kidney,kidney.pred)
proc.kidney$auc
ggroc(proc.kidney)


#### other pick of variable
library(hash)
load("filter_data.RData")
###### bone
bone.df <- cbind(data.frame(inter_candidate$bone),bone = as.factor(bone))
###### brain
brain.df <- cbind(data.frame(inter_candidate$brain),brain = as.factor(brain))
###### lung
lung.df <- cbind(data.frame(inter_candidate$lung),lung = as.factor(lung))
##### liver
liver.df <- cbind(data.frame(inter_candidate$liver),liver = as.factor(liver))
##### kidney
kidney.df <- cbind(data.frame(inter_candidate$kidney),kidney = as.factor(kidney))


#### by 器官
library(pROC)
#### bone
set.seed(2000)
train.bone <- dplyr::slice_sample(bone.df,prop = 0.8,by = bone)
test.bone <- dplyr::anti_join(bone.df,train.bone)
table(train.bone$bone)
table(test.bone$bone)
logit.train.bone <- glm(bone~.,data = train.bone,family = 'binomial')
pred.test.bone <- predict(logit.train.bone,type = 'response',newdata = test.bone )
pred.bone <- ifelse(pred.test.bone>0.5,1,0)
bone_table <- table(predicted = pred.bone,actual = test.bone$bone)
bone_table
bone.pred <- predict(logit.train.bone,type = "response")
proc.bone <- roc(train.bone$bone,bone.pred)
proc.bone$auc
ggroc(proc.bone)



logit.bone <- glm(bone~.,data = bone.df,family = "binomial")
summary(logit.bone)
bone.pred <- predict(logit.bone,type = "response")
proc.bone <- roc(bone,bone.pred)
proc.bone$auc
ggroc(proc.bone)


##### brain
logit.brain <- glm(brain~.,data = brain.df,family = "binomial")
summary(logit.brain)
brain.pred <- predict(logit.brain,type = "response")
proc.brain <- roc(brain,brain.pred)
proc.brain$auc
ggroc(proc.brain)

#### lung
logit.lung <- glm(lung~.,data = lung.df,family = "binomial")
summary(logit.lung)
lung.pred <- predict(logit.lung,type = "response")
proc.lung <- roc(lung,lung.pred)
proc.lung$auc
ggroc(proc.lung)

#### liver
logit.liver <- glm(liver~.,data = liver.df,family = "binomial")
summary(logit.liver)
liver.pred <- predict(logit.liver,type = "response")
proc.liver <- roc(liver,liver.pred)
proc.liver$auc
ggroc(proc.liver)


#### kidney
logit.kidney <- glm(kidney~.,data = kidney.df,family = "binomial")
summary(logit.kidney)
kidney.pred <- predict(logit.kidney,type = "response")
proc.kidney <- roc(kidney,kidney.pred)
proc.kidney$auc
ggroc(proc.kidney)

