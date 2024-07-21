## ----Read Data
library(tidyverse)
library(stargazer)
library(pscl)
library(pROC)
library(hash)
load("glm_model/Datas.RData")
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
## ----Count bar plot
barplot(table(count_response_label),main = "Count response Barplot")
## ----create Interaction
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

## ----Create Interaction DataFrame(lassohurdle)
load("glm_model/binary.Rdata")
load("glm_model/count.Rdata")
load("glm_model/response.Rdata")
load("glm_model/variable.Rdata")
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
## ----Zero Inflated(lasso)
hurdlemodel <- hurdle(as.formula(paste("count_response_label~",paste(selected_variables_count[2:length(selected_variables_count)],collapse  = "+"),"|",paste(selected_variables_binary[2:length(selected_variables_binary)],collapse = "+"))),data = df.analysis)
summary(hurdlemodel)
## ----Zero Inflated(lasso) consequence
p <- ifelse(round(predict(hurdlemodel))>5,5,round(predict(hurdlemodel)))
MSE <- mean((p-count_response_label)^2)
var.count <- var(count_response_label)
knitr::kable(data.frame(MSE = c(MSE),variance = c(var.count)))

## ----Create Interaction DataFrame(lassologit)
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
## ----Logistic(lasso)
#### bone
logit.bone <- glm(bone~.,data = bone.df,family = "binomial")
#summary(logit.bone)
bone.pred <- predict(logit.bone,type = "response")
proc.bone <- roc(bone,bone.pred)
#proc.bone$auc
#ggroc(proc.bone)


##### brain
logit.brain <- glm(brain~.,data = brain.df,family = "binomial")
#summary(logit.brain)
brain.pred <- predict(logit.brain,type = "response")
proc.brain <- roc(brain,brain.pred)
#proc.brain$auc
#ggroc(proc.brain)

#### lung
logit.lung <- glm(lung~.,data = lung.df,family = "binomial",maxit=2000)
#summary(logit.lung)
lung.pred <- predict(logit.lung,type = "response")
proc.lung <- roc(lung,lung.pred)
#proc.lung$auc
#ggroc(proc.lung)

#### liver
logit.liver <- glm(liver~.,data = liver.df,family = "binomial",maxit=2000)
#summary(logit.liver)
liver.pred <- predict(logit.liver,type = "response")
proc.liver <- roc(liver,liver.pred)
#proc.liver$auc
#ggroc(proc.liver)


#### kidney
logit.kidney <- glm(kidney~.,data = kidney.df,family = "binomial")
#summary(logit.kidney)
kidney.pred <- predict(logit.kidney,type = "response")
proc.kidney <- roc(kidney,kidney.pred)
#proc.kidney$auc
#ggroc(proc.kidney)
## ----overall summary(logit.lasso)
stargazer(logit.bone,logit.brain,logit.kidney,logit.liver,logit.lung,type = 'html')
## ----auc(logit.lasso)
auc.df <- data.frame(bone = c(proc.bone$auc),
                     brain = c(proc.brain$auc),
                     liver = c(proc.liver$auc),
                     lung = c(proc.lung$auc),
                     kidney = c(proc.kidney$auc))
knitr::kable(auc.df)

## ----Train Test
## bone
set.seed(2000)
train.bone <- dplyr::slice_sample(bone.df,prop = 0.8,by = bone)
test.bone <- dplyr::anti_join(bone.df,train.bone)
logit.train.bone <- glm(bone~.,data = train.bone,family = 'binomial')
pred.test.bone <- predict(logit.train.bone,type = 'response',newdata = test.bone )
pred.bone <- ifelse(pred.test.bone>0.5,1,0)
bone_table <- table(predicted = pred.bone,actual = test.bone$bone)
bone_acc <- (bone_table[1,1]+bone_table[2,2])/sum(bone_table)

## brain
set.seed(2000)
train.brain <- dplyr::slice_sample(brain.df,prop = 0.8,by = brain)
test.brain <- dplyr::anti_join(brain.df,train.brain)
logit.train.brain <- glm(brain~.,data = train.brain,family = 'binomial')
pred.test.brain <- predict(logit.train.brain,type = 'response',newdata = test.brain )
pred.brain <- ifelse(pred.test.brain>0.5,1,0)
brain_table <- table(predicted = pred.brain,actual = test.brain$brain)
brain_acc <- (brain_table[1,1]+brain_table[2,2])/sum(brain_table)


## liver
set.seed(2000)
train.liver <- dplyr::slice_sample(liver.df,prop = 0.8,by = liver)
test.liver <- dplyr::anti_join(liver.df,train.liver)
logit.train.liver <- glm(liver~.,data = train.liver,family = 'binomial',maxit=1000)
pred.test.liver <- predict(logit.train.liver,type = 'response',newdata = test.liver )
pred.liver <- ifelse(pred.test.liver>0.5,1,0)
liver_table <- table(predicted = pred.liver,actual = test.liver$liver)
liver_acc <- (liver_table[1,1]+liver_table[2,2])/sum(liver_table)


## lung
set.seed(2000)
train.lung <- dplyr::slice_sample(lung.df,prop = 0.8,by = lung)
test.lung <- dplyr::anti_join(lung.df,train.lung)
logit.train.lung <- glm(lung~.,data = train.lung,family = 'binomial',maxit = 1000)
pred.test.lung <- predict(logit.train.lung,type = 'response',newdata = test.lung )
pred.lung <- ifelse(pred.test.lung>0.5,1,0)
lung_table <- table(predicted = pred.lung,actual = test.lung$lung)
lung_acc <- (lung_table[1,1]+lung_table[2,2])/sum(lung_table)


## kidney
set.seed(2000)
train.kidney <- dplyr::slice_sample(kidney.df,prop = 0.8,by = kidney)
test.kidney <- dplyr::anti_join(kidney.df,train.kidney)
logit.train.kidney <- glm(kidney~.,data = train.kidney,family = 'binomial')
pred.test.kidney <- predict(logit.train.kidney,type = 'response',newdata = test.kidney )
pred.kidney <- ifelse(pred.test.kidney>0.5,1,0)
kidney_table <- table(predicted = pred.kidney,actual = test.kidney$kidney)
kidney_acc <- (kidney_table[1,1]+kidney_table[2,2])/sum(kidney_table)

acc.df <- data.frame(bone = c(bone_acc),
                     brain = c(brain_acc),
                     liver = c(liver_acc),
                     lung = c(lung_acc),
                     kidney = c(kidney_acc))
knitr::kable(acc.df)

## ----Zero Inflated(pvalue)
load("glm_model/filter_data2.RData")
organ_candidate <- data.frame(organ_candidate) %>%
  mutate(binary_response_label,count_response_label)

selected_variables <- sapply(colnames(organ_candidate)[1:(length(colnames(organ_candidate))-2)],add_symbol)
hurdlemodel <- hurdle(as.formula(paste("count_response_label~",paste(selected_variables,collapse  = "+"),"|",paste(selected_variables,collapse = "+"))),data = organ_candidate)
summary(hurdlemodel)
## ----Zero Inflated(pvalue) consequence
p <- ifelse(round(predict(hurdlemodel))>5,5,round(predict(hurdlemodel)))
MSE <- mean((p-count_response_label)^2)
var.count <- var(count_response_label)
knitr::kable(data.frame(MSE = c(MSE),variance = c(var.count)))
## ----Logistic(pvalue)
load("glm_model/filter_data.RData")
library(stringr)
remove_special_chars <- function(covar_labels){
  covar_labels %>% 
    str_replace_all("_",   ".")
}

###### bone
bone.df <- cbind(data.frame(inter_candidate$bone),bone = as.factor(bone))
colnames(bone.df) <- remove_special_chars(colnames(bone.df))
###### brain
brain.df <- cbind(data.frame(inter_candidate$brain),brain = as.factor(brain))
colnames(brain.df) <- remove_special_chars(colnames(brain.df))
###### lung
lung.df <- cbind(data.frame(inter_candidate$lung),lung = as.factor(lung))
colnames(lung.df) <- remove_special_chars(colnames(lung.df))
##### liver
liver.df <- cbind(data.frame(inter_candidate$liver),liver = as.factor(liver))
colnames(liver.df) <- remove_special_chars(colnames(liver.df))
##### kidney
kidney.df <- cbind(data.frame(inter_candidate$kidney),kidney = as.factor(kidney))
colnames(kidney.df) <- remove_special_chars(colnames(kidney.df))

#### by 器官
#### bone
logit.bone <- glm(bone~.,data = bone.df,family = "binomial")
#summary(logit.bone)
bone.pred <- predict(logit.bone,type = "response")
proc.bone <- roc(bone,bone.pred)
#proc.bone$auc
#ggroc(proc.bone)


##### brain
logit.brain <- glm(brain~.,data = brain.df,family = "binomial")
#summary(logit.brain)
brain.pred <- predict(logit.brain,type = "response")
proc.brain <- roc(brain,brain.pred)
#proc.brain$auc
#ggroc(proc.brain)

#### lung
logit.lung <- glm(lung~.,data = lung.df,family = "binomial")
#summary(logit.lung)
lung.pred <- predict(logit.lung,type = "response")
proc.lung <- roc(lung,lung.pred)
#proc.lung$auc
#ggroc(proc.lung)

#### liver
logit.liver <- glm(liver~.,data = liver.df,family = "binomial")
#summary(logit.liver)
liver.pred <- predict(logit.liver,type = "response")
proc.liver <- roc(liver,liver.pred)
#proc.liver$auc
#ggroc(proc.liver)


#### kidney
logit.kidney <- glm(kidney~.,data = kidney.df,family = "binomial")
#summary(logit.kidney)
kidney.pred <- predict(logit.kidney,type = "response")
proc.kidney <- roc(kidney,kidney.pred)
#proc.kidney$auc
#ggroc(proc.kidney)
## ----overall summary(logit.pv)
stargazer(logit.bone,logit.brain,logit.kidney,logit.liver,logit.lung,type = 'html')
## ----auc(logit.pv)
auc.df <- data.frame(bone = c(proc.bone$auc),
                     brain = c(proc.brain$auc),
                     liver = c(proc.liver$auc),
                     lung = c(proc.lung$auc),
                     kidney = c(proc.kidney$auc))
knitr::kable(auc.df)

## ----Train Test(p-value)
## bone
set.seed(2000)
train.bone <- dplyr::slice_sample(bone.df,prop = 0.8,by = bone)
test.bone <- dplyr::anti_join(bone.df,train.bone)
logit.train.bone <- glm(bone~.,data = train.bone,family = 'binomial')
pred.test.bone <- predict(logit.train.bone,type = 'response',newdata = test.bone )
pred.bone <- ifelse(pred.test.bone>0.5,1,0)
bone_table <- table(predicted = pred.bone,actual = test.bone$bone)
bone_acc <- (bone_table[1,1]+bone_table[2,2])/sum(bone_table)

## brain
set.seed(2000)
train.brain <- dplyr::slice_sample(brain.df,prop = 0.8,by = brain)
test.brain <- dplyr::anti_join(brain.df,train.brain)
logit.train.brain <- glm(brain~.,data = train.brain,family = 'binomial')
pred.test.brain <- predict(logit.train.brain,type = 'response',newdata = test.brain )
pred.brain <- ifelse(pred.test.brain>0.5,1,0)
brain_table <- table(predicted = pred.brain,actual = test.brain$brain)
brain_acc <- (brain_table[1,1]+brain_table[2,2])/sum(brain_table)


## liver
set.seed(2000)
train.liver <- dplyr::slice_sample(liver.df,prop = 0.8,by = liver)
test.liver <- dplyr::anti_join(liver.df,train.liver)
logit.train.liver <- glm(liver~.,data = train.liver,family = 'binomial',maxit=1000)
pred.test.liver <- predict(logit.train.liver,type = 'response',newdata = test.liver )
pred.liver <- ifelse(pred.test.liver>0.5,1,0)
liver_table <- table(predicted = pred.liver,actual = test.liver$liver)
liver_acc <- (liver_table[1,1]+liver_table[2,2])/sum(liver_table)


## lung
set.seed(2000)
train.lung <- dplyr::slice_sample(lung.df,prop = 0.8,by = lung)
test.lung <- dplyr::anti_join(lung.df,train.lung)
logit.train.lung <- glm(lung~.,data = train.lung,family = 'binomial',maxit = 1000)
pred.test.lung <- predict(logit.train.lung,type = 'response',newdata = test.lung )
pred.lung <- ifelse(pred.test.lung>0.5,1,0)
lung_table <- table(predicted = pred.lung,actual = test.lung$lung)
lung_acc <- (lung_table[1,1]+lung_table[2,2])/sum(lung_table)


## kidney
set.seed(2000)
train.kidney <- dplyr::slice_sample(kidney.df,prop = 0.8,by = kidney)
test.kidney <- dplyr::anti_join(kidney.df,train.kidney)
logit.train.kidney <- glm(kidney~.,data = train.kidney,family = 'binomial')
pred.test.kidney <- predict(logit.train.kidney,type = 'response',newdata = test.kidney )
pred.kidney <- ifelse(pred.test.kidney>0.5,1,0)
kidney_table <- table(predicted = pred.kidney,actual = test.kidney$kidney)
kidney_acc <- (kidney_table[1,1]+kidney_table[2,2])/sum(kidney_table)

acc.df <- data.frame(bone = c(bone_acc),
                     brain = c(brain_acc),
                     liver = c(liver_acc),
                     lung = c(lung_acc),
                     kidney = c(kidney_acc))
knitr::kable(acc.df)