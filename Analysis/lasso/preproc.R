# 匯入資料
data_genes = read.csv("C:/Users/user/Downloads/genes.csv")
data_meta = read.csv("C:/Users/user/Documents/成功大學/統計諮詢/group3_final_report/Analysis/Data/metabolites.csv")
dim(data_genes) # 1406 19222
dim(data_meta) # 265 226

# outer join
# merged_data <- merge(data_genes, data_meta, by.x  = "X", by.y  = "ID", all = TRUE)
# merged_data <- merged_data[, !(names(merged_data) %in% c("ID", "X"))]
# dim(merged_data) # 1447 19446

# inner join (genes, metabolites)
merged_data2 <- merge(data_genes, data_meta, by.x  = "X", by.y  = "ID", all = FALSE)
dim(merged_data2) # 265 19447
merged_data2 <- merged_data2[, !(names(merged_data2) %in% c("X", "ID"))]
dim(merged_data2) # 265 19446
str(merged_data2)



# check 不為 num 的變數
var_types <- sapply(merged_data2, class)
non_numeric_vars <- var_types[var_types != "numeric"]
non_numeric_vars # 所有變數都是num

# check NA
na_counts_row <- apply(merged_data2, 1, function(x) sum(is.na(x)))
na_counts_row # no NA

# 移除零變異行 對pca沒有貢獻
zero = which(apply(merged_data2, 2, sd) == 0)
merged_data2_cleaned <- merged_data2[, -zero]
dim(merged_data2_cleaned)

# PCA
pca_gm <- prcomp(merged_data2_cleaned, center = TRUE, scale=TRUE)
pve_gm = (pca_gm$sdev^2)/sum(pca_gm$sdev^2) #* 100
summary(pca_gm)
# Plot
plot(pve_gm, type="o", ylab="Eigenvalue", xlab = "Principal Component", main='genes + metabolites')
plot(cumsum(pve_gm), type="o", ylab="Explained Variance ", xlab = "Principal Component", main='genes + metabolites')



# -------------只拿 gene 資料 PCA ----------------------
dim(data_genes) # 1406 19222
na_counts_row <- apply(data_genes, 1, function(x) sum(is.na(x)))
which(na_counts_row!=0) # 每一列皆沒 NA

# 移除ID
data_genes = data_genes[, !names(data_genes) %in% "X"]

# PCA
pca_g <- prcomp(data_genes, center = TRUE, scale=TRUE)
pve_g = (pca_g$sdev^2)/sum(pca_g$sdev^2) #* 100

# Plot
plot(pve_g , type="o", ylab="Eigenvalue", xlab = "Principal Component", main='genes')
plot(cumsum(pve_g), type="o", ylab="Explained Variance ", xlab = "Principal Component", main='genes')

#-------------只拿metabolites做pca-----------------
data_meta = data_meta[, !names(data_meta) %in% "ID"]
pca_m <- prcomp(data_meta, center = TRUE, scale=TRUE)
pve_m = (pca_m$sdev^2)/sum(pca_m$sdev^2) #* 100

plot(pve_m , type="o", ylab="Eigenvalue", xlab = "Principal Component", main='metabolites')
plot(cumsum(pve_m), type="o", ylab="Explained Variance ", xlab = "Principal Component", main='metabolites')



##-------------Lasso(meta+gene)-----------
data_response = read.csv("C:/Users/user/Documents/成功大學/統計諮詢/group3_final_report/Analysis/Data/Response.csv")

data_y <- data_response[match(merged_data2$X, data_response$ID), ]

#bone
data_y$bone <- ifelse(data_y$bone == "Metastasis", 1, 0)

install.packages("glmnet")
library(glmnet)

merged_data2_scale <- scale(merged_data2_cleaned)

lasso_logistic_bone <- glmnet(merged_data2_scale,data_y$bone, family = "binomial", alpha = 1)
cv_lasso_bone <- cv.glmnet(merged_data2_scale,data_y$bone, family = "binomial", alpha = 1)

final_lasso_bone <- glmnet(merged_data2_scale, data_y$bone, family = "binomial", alpha = 1, lambda = cv_lasso_bone$lambda.min)

lasso_coefficients_bone <- coef(final_lasso_bone)
coefficients_matrix_bone <- as.matrix(lasso_coefficients_bone)
selected_variables_bone <- rownames(coefficients_matrix_bone)[coefficients_matrix_bone != 0]
selected_variables_bone

#brain
data_y$brain <- ifelse(data_y$brain == "Metastasis", 1, 0)  

lasso_logistic_brain <- glmnet(merged_data2_scale,data_y$brain, family = "binomial", alpha = 1)
cv_lasso_brain <- cv.glmnet(merged_data2_scale,data_y$brain, family = "binomial", alpha = 1)

final_lasso_brain <- glmnet(merged_data2_scale, data_y$brain, family = "binomial", alpha = 1, lambda = cv_lasso_brain$lambda.min)

lasso_coefficients_brain <- coef(final_lasso_brain)
coefficients_matrix_brain <- as.matrix(lasso_coefficients_brain)
selected_variables_brain <- rownames(coefficients_matrix_brain)[coefficients_matrix_brain != 0]
selected_variables_brain

#kidney
data_y$kidney <- ifelse(data_y$kidney == "Metastasis", 1, 0)

lasso_logistic_kidney <- glmnet(merged_data2_scale,data_y$kidney, family = "binomial", alpha = 1)
cv_lasso_kidney <- cv.glmnet(merged_data2_scale,data_y$kidney, family = "binomial", alpha = 1)

final_lasso_kidney <- glmnet(merged_data2_scale, data_y$kidney, family = "binomial", alpha = 1, lambda = cv_lasso_kidney$lambda.min)

lasso_coefficients_kidney <- coef(final_lasso_kidney)
coefficients_matrix_kidney <- as.matrix(lasso_coefficients_kidney)
selected_variables_kidney <- rownames(coefficients_matrix_kidney)[coefficients_matrix_kidney != 0]
selected_variables_kidney

#liver
data_y$liver <- ifelse(data_y$liver == "Metastasis", 1, 0)

lasso_logistic_liver <- glmnet(merged_data2_scale,data_y$liver, family = "binomial", alpha = 1)
cv_lasso_liver <- cv.glmnet(merged_data2_scale,data_y$liver, family = "binomial", alpha = 1)

final_lasso_liver <- glmnet(merged_data2_scale, data_y$liver, family = "binomial", alpha = 1, lambda = cv_lasso_liver$lambda.min)

lasso_coefficients_liver <- coef(final_lasso_liver )
coefficients_matrix_liver  <- as.matrix(lasso_coefficients_liver )
selected_variables_liver  <- rownames(coefficients_matrix_liver )[coefficients_matrix_liver != 0]
selected_variables_liver 

#lung
data_y$lung <- ifelse(data_y$lung == "Metastasis", 1, 0)

lasso_logistic_lung <- glmnet(merged_data2_scale,data_y$lung, family = "binomial", alpha = 1)
cv_lasso_lung <- cv.glmnet(merged_data2_scale,data_y$lung, family = "binomial", alpha = 1)

final_lasso_lung <- glmnet(merged_data2_scale, data_y$lung, family = "binomial", alpha = 1, lambda = cv_lasso_lung$lambda.min)

lasso_coefficients_lung <- coef(final_lasso_lung )
coefficients_matrix_lung  <- as.matrix(lasso_coefficients_lung)
selected_variables_lung  <- rownames(coefficients_matrix_lung)[coefficients_matrix_lung != 0]
selected_variables_lung



##基因加代謝體對response做lasso(只要一個轉移就算轉移)
which(table(data_y$ID)>1)
data_y <- unique(data_y)
merged_data2_scale <- unique(merged_data2_scale)

data_y$sum <- data_y$bone+data_y$brain+data_y$kidney+data_y$liver+data_y$lung
data_y
data_y$sum <- ifelse(data_y$sum > 0, 1, 0)

lasso_logistic_sum <- glmnet(merged_data2_scale,data_y$sum, family = "binomial", alpha = 1)
cv_lasso_sum <- cv.glmnet(merged_data2_scale,data_y$sum, family = "binomial", alpha = 1)

cv <- c()
for (i in 1:100){
  set.seed(i*100)
  cv_lasso_sum <- cv.glmnet(merged_data2_scale,data_y$sum, family = "binomial", alpha = 1)
  cv=append(cv,cv_lasso_sum$lambda.min)
  print(i)
}
cv

final_lasso_sum <- glmnet(merged_data2_scale, data_y$sum, family = "binomial", alpha = 1, lambda = min(cv))

lasso_coefficients_sum <- coef(final_lasso_sum)
coefficients_matrix_sum  <- as.matrix(lasso_coefficients_sum )
selected_variables_sum  <- rownames(coefficients_matrix_sum)[coefficients_matrix_sum != 0]
selected_variables_sum

save(selected_variables_sum, file = "C:/Users/user/Desktop/sum.Rdata")


#前處理
load("Datas.RData")

#######對每個癌症各挑25個變數
load("interaction_transfer.Rdata")

Response_filter$bone <- ifelse(Response_filter$bone == "Metastasis", 1, 0)
Response_filter$brain <- ifelse(Response_filter$brain == "Metastasis", 1, 0)
Response_filter$kidney <- ifelse(Response_filter$kidney == "Metastasis", 1, 0)
Response_filter$liver <- ifelse(Response_filter$liver == "Metastasis", 1, 0)
Response_filter$lung <- ifelse(Response_filter$lung == "Metastasis", 1, 0)

install.packages("glmnet")
library(glmnet)

genes_filter_scale <- scale(genes_filter)
metabolites_filter_scale <- scale(metabolites_filter)


##bone

#基因
final_genes_bone <- glmnet(genes_filter_scale, Response_filter$bone, family = "binomial", alpha = 1, lambda = 0.076)
coefficients_genes_bone <- coef(final_genes_bone)
coefficients_matrix_genes_bone <- as.matrix(coefficients_genes_bone)
variables_genes_bone <- rownames(coefficients_matrix_genes_bone )[coefficients_matrix_genes_bone  != 0]
variables_genes_bone

#代謝體
final_metabolites_bone <- glmnet(metabolites_filter_scale, Response_filter$bone, family = "binomial", alpha = 1, lambda = 0.047)
coefficients_metabolites_bone <- coef(final_metabolites_bone)
coefficients_matrix_metabolites_bone <- as.matrix(coefficients_metabolites_bone)
variables_metabolites_bone <- rownames(coefficients_matrix_metabolites_bone )[coefficients_matrix_metabolites_bone  != 0]
variables_metabolites_bone

#考慮交互作用再用lasso
variables_genes_bone <- variables_genes_bone[variables_genes_bone != "(Intercept)"]
genes_bone.df <- genes_filter[, variables_genes_bone]

variables_metabolites_bone <- variables_metabolites_bone[variables_metabolites_bone != "(Intercept)"]
metabolites_bone.df <- metabolites_filter[, variables_metabolites_bone]

x_bone <- cbind(genes_bone.df , metabolites_bone.df)

bone_interaction <- interaction.trasfer.function(x_bone)
bone_interaction_scale <- scale(bone_interaction)
lasso_bone <- glmnet(bone_interaction_scale, Response_filter$bone, family = "binomial", alpha = 1, lambda = 0.02)
coefficients_lasso_bone <- coef(lasso_bone)
coefficients_matrix_lasso_bone <- as.matrix(coefficients_lasso_bone)
variables_bone <- rownames(coefficients_matrix_lasso_bone )[coefficients_matrix_lasso_bone  != 0]
variables_bone


##brain

#基因
final_genes_brain <- glmnet(genes_filter_scale, Response_filter$brain, family = "binomial", alpha = 1, lambda = 0.07)
coefficients_genes_brain <- coef(final_genes_brain)
coefficients_matrix_genes_brain <- as.matrix(coefficients_genes_brain)
variables_genes_brain <- rownames(coefficients_matrix_genes_brain )[coefficients_matrix_genes_brain  != 0]
variables_genes_brain

#代謝體
final_metabolites_brain <- glmnet(metabolites_filter_scale, Response_filter$brain, family = "binomial", alpha = 1, lambda = 0.04)
coefficients_metabolites_brain <- coef(final_metabolites_brain)
coefficients_matrix_metabolites_brain <- as.matrix(coefficients_metabolites_brain)
variables_metabolites_brain <- rownames(coefficients_matrix_metabolites_brain )[coefficients_matrix_metabolites_brain  != 0]
variables_metabolites_brain

#考慮交互作用再用lasso
variables_genes_brain <- variables_genes_brain[variables_genes_brain != "(Intercept)"]
genes_brain.df <- genes_filter[, variables_genes_brain]

variables_metabolites_brain <- variables_metabolites_brain[variables_metabolites_brain != "(Intercept)"]
metabolites_brain.df <- metabolites_filter[, variables_metabolites_brain]

x_brain <- cbind(genes_brain.df , metabolites_brain.df)

brain_interaction <- interaction.trasfer.function(x_brain)
brain_interaction_scale <- scale(brain_interaction)
lasso_brain <- glmnet(brain_interaction_scale, Response_filter$brain, family = "binomial", alpha = 1, lambda = 0.02)
coefficients_lasso_brain <- coef(lasso_brain)
coefficients_matrix_lasso_brain <- as.matrix(coefficients_lasso_brain)
variables_brain <- rownames(coefficients_matrix_lasso_brain )[coefficients_matrix_lasso_brain  != 0]
variables_brain

##kidney

#基因
final_genes_kidney <- glmnet(genes_filter_scale, Response_filter$kidney, family = "binomial", alpha = 1, lambda = 0.073)
coefficients_genes_kidney <- coef(final_genes_kidney)
coefficients_matrix_genes_kidney <- as.matrix(coefficients_genes_kidney)
variables_genes_kidney <- rownames(coefficients_matrix_genes_kidney )[coefficients_matrix_genes_kidney  != 0]
variables_genes_kidney

#代謝體
final_metabolites_kidney <- glmnet(metabolites_filter_scale, Response_filter$kidney, family = "binomial", alpha = 1, lambda = 0.04)
coefficients_metabolites_kidney <- coef(final_metabolites_kidney)
coefficients_matrix_metabolites_kidney <- as.matrix(coefficients_metabolites_kidney)
variables_metabolites_kidney <- rownames(coefficients_matrix_metabolites_kidney )[coefficients_matrix_metabolites_kidney  != 0]
variables_metabolites_kidney

#考慮交互作用再用lasso
variables_genes_kidney <- variables_genes_kidney[variables_genes_kidney != "(Intercept)"]
genes_kidney.df <- genes_filter[, variables_genes_kidney]

variables_metabolites_kidney <- variables_metabolites_kidney[variables_metabolites_kidney != "(Intercept)"]
metabolites_kidney.df <- metabolites_filter[, variables_metabolites_kidney]

x_kidney <- cbind(genes_kidney.df , metabolites_kidney.df)

kidney_interaction <- interaction.trasfer.function(x_kidney)
kidney_interaction_scale <- scale(kidney_interaction)
lasso_kidney <- glmnet(kidney_interaction_scale, Response_filter$kidney, family = "binomial", alpha = 1, lambda = 0.014)
coefficients_lasso_kidney <- coef(lasso_kidney)
coefficients_matrix_lasso_kidney <- as.matrix(coefficients_lasso_kidney)
variables_kidney <- rownames(coefficients_matrix_lasso_kidney )[coefficients_matrix_lasso_kidney  != 0]
variables_kidney

##liver

#基因
final_genes_liver <- glmnet(genes_filter_scale, Response_filter$liver, family = "binomial", alpha = 1, lambda = 0.073)
coefficients_genes_liver <- coef(final_genes_liver)
coefficients_matrix_genes_liver <- as.matrix(coefficients_genes_liver)
variables_genes_liver <- rownames(coefficients_matrix_genes_liver )[coefficients_matrix_genes_liver  != 0]
variables_genes_liver

#代謝體
final_metabolites_liver <- glmnet(metabolites_filter_scale, Response_filter$liver, family = "binomial", alpha = 1, lambda = 0.022)
coefficients_metabolites_liver <- coef(final_metabolites_liver)
coefficients_matrix_metabolites_liver <- as.matrix(coefficients_metabolites_liver)
variables_metabolites_liver <- rownames(coefficients_matrix_metabolites_liver )[coefficients_matrix_metabolites_liver  != 0]
variables_metabolites_liver

#考慮交互作用再用lasso
variables_genes_liver <- variables_genes_liver[variables_genes_liver != "(Intercept)"]
genes_liver.df <- genes_filter[, variables_genes_liver]

variables_metabolites_liver <- variables_metabolites_liver[variables_metabolites_liver != "(Intercept)"]
metabolites_liver.df <- metabolites_filter[, variables_metabolites_liver]

x_liver <- cbind(genes_liver.df , metabolites_liver.df)

liver_interaction <- interaction.trasfer.function(x_liver)
liver_interaction_scale <- scale(liver_interaction)
lasso_liver <- glmnet(liver_interaction_scale, Response_filter$liver, family = "binomial", alpha = 1, lambda = 0.028)
coefficients_lasso_liver <- coef(lasso_liver)
coefficients_matrix_lasso_liver <- as.matrix(coefficients_lasso_liver)
variables_liver <- rownames(coefficients_matrix_lasso_liver )[coefficients_matrix_lasso_liver != 0]
variables_liver

##lung

#基因
final_genes_lung <- glmnet(genes_filter_scale, Response_filter$lung, family = "binomial", alpha = 1, lambda = 0.052)
coefficients_genes_lung <- coef(final_genes_lung)
coefficients_matrix_genes_lung <- as.matrix(coefficients_genes_lung)
variables_genes_lung <- rownames(coefficients_matrix_genes_lung )[coefficients_matrix_genes_lung  != 0]
variables_genes_lung

#代謝體
final_metabolites_lung <- glmnet(metabolites_filter_scale, Response_filter$lung, family = "binomial", alpha = 1, lambda = 0.022)
coefficients_metabolites_lung <- coef(final_metabolites_lung)
coefficients_matrix_metabolites_lung <- as.matrix(coefficients_metabolites_lung)
variables_metabolites_lung <- rownames(coefficients_matrix_metabolites_lung )[coefficients_matrix_metabolites_lung  != 0]
variables_metabolites_lung

#考慮交互作用再用lasso
variables_genes_lung <- variables_genes_lung[variables_genes_lung != "(Intercept)"]
genes_lung.df <- genes_filter[, variables_genes_lung]

variables_metabolites_lung <- variables_metabolites_lung[variables_metabolites_lung != "(Intercept)"]
metabolites_lung.df <- metabolites_filter[, variables_metabolites_lung]

x_lung <- cbind(genes_lung.df , metabolites_lung.df)

lung_interaction <- interaction.trasfer.function(x_lung)
lung_interaction_scale <- scale(lung_interaction)
lasso_lung <- glmnet(lung_interaction_scale, Response_filter$lung, family = "binomial", alpha = 1, lambda = 0.07)
coefficients_lasso_lung <- coef(lasso_lung)
coefficients_matrix_lasso_lung <- as.matrix(coefficients_lasso_lung)
variables_lung <- rownames(coefficients_matrix_lasso_lung)[coefficients_matrix_lasso_lung != 0]
variables_lung

##存檔
save(variables_bone, variables_brain, variables_kidney, variables_liver, variables_lung, file = "C:/Users/user/Desktop/response.Rdata")


######有一個轉移就算轉移
Response_filter$sum <- Response_filter$bone+Response_filter$brain+Response_filter$kidney+Response_filter$liver+Response_filter$lung
Response_filter
Response_filter$sum <- ifelse(Response_filter$sum > 0, 1, 0)

#基因
final_genes_sum <- glmnet(genes_filter_scale, Response_filter$sum, family = "binomial", alpha = 1, lambda = 0.044)
coefficients_genes_sum <- coef(final_genes_sum)
coefficients_matrix_genes_sum <- as.matrix(coefficients_genes_sum)
variables_genes_sum <- rownames(coefficients_matrix_genes_sum )[coefficients_matrix_genes_sum  != 0]
variables_genes_sum

#代謝體
final_metabolites_sum <- glmnet(metabolites_filter_scale, Response_filter$sum, family = "binomial", alpha = 1, lambda = 0.02)
coefficients_metabolites_sum <- coef(final_metabolites_sum)
coefficients_matrix_metabolites_sum <- as.matrix(coefficients_metabolites_sum)
variables_metabolites_sum <- rownames(coefficients_matrix_metabolites_sum )[coefficients_matrix_metabolites_sum != 0]
variables_metabolites_sum

#考慮交互作用再用lasso
variables_genes_sum <- variables_genes_sum[variables_genes_sum != "(Intercept)"]
genes_sum.df <- genes_filter[, variables_genes_sum]

variables_metabolites_sum <- variables_metabolites_sum[variables_metabolites_sum != "(Intercept)"]
metabolites_sum.df <- metabolites_filter[, variables_metabolites_sum]

x_sum <- cbind(genes_sum.df , metabolites_sum.df)

sum_interaction <- interaction.trasfer.function(x_sum)
sum_interaction_scale <- scale(sum_interaction)
lasso_sum <- glmnet(sum_interaction_scale, Response_filter$sum, family = "binomial", alpha = 1, lambda = 0.05)
coefficients_lasso_sum <- coef(lasso_sum)
coefficients_matrix_lasso_sum <- as.matrix(coefficients_lasso_sum)
variables_binary <- rownames(coefficients_matrix_lasso_sum)[coefficients_matrix_lasso_sum != 0]
variables_binary

#存檔
save(variables_binary, file = "C:/Users/user/Desktop/binary.Rdata")


#######轉移幾個地方
Response_filter$count <- Response_filter$bone+Response_filter$brain+Response_filter$kidney+Response_filter$liver+Response_filter$lung

#基因
final_genes_count <- glmnet(genes_filter_scale, Response_filter$count, family = "poisson", alpha = 1, lambda = 0.21)
coefficients_genes_count <- coef(final_genes_count)
coefficients_matrix_genes_count <- as.matrix(coefficients_genes_count)
variables_genes_count <- rownames(coefficients_matrix_genes_count )[coefficients_matrix_genes_count  != 0]
variables_genes_count

#代謝體
final_metabolites_count <- glmnet(metabolites_filter_scale, Response_filter$count, family = "poisson", alpha = 1, lambda = 0.092)
coefficients_metabolites_count <- coef(final_metabolites_count)
coefficients_matrix_metabolites_count <- as.matrix(coefficients_metabolites_count)
variables_metabolites_count <- rownames(coefficients_matrix_metabolites_count)[coefficients_matrix_metabolites_count != 0]
variables_metabolites_count

#考慮交互作用再用lasso
variables_genes_count <- variables_genes_count[variables_genes_count != "(Intercept)"]
genes_count.df <- genes_filter[, variables_genes_count]

variables_metabolites_count <- variables_metabolites_count[variables_metabolites_count != "(Intercept)"]
metabolites_count.df <- metabolites_filter[, variables_metabolites_count]

x_count <- cbind(genes_count.df , metabolites_count.df)

count_interaction <- interaction.trasfer.function(x_count)
count_interaction_scale <- scale(count_interaction)
lasso_count <- glmnet(count_interaction_scale, Response_filter$count, family = "poisson", alpha = 1, lambda = 0.1)
coefficients_lasso_count <- coef(lasso_count)
coefficients_matrix_lasso_count <- as.matrix(coefficients_lasso_count)
variables_count <- rownames(coefficients_matrix_lasso_count)[coefficients_matrix_lasso_count != 0]
variables_count

#存檔
save(variables_count, file = "C:/Users/user/Desktop/count.Rdata")

which(variables_count== "ABCA7..10347. _ CD44..960." )


x_binary <- x_sum
save(x_bone, x_brain, x_kidney, x_liver, x_lung, x_binary, x_count, file = "C:/Users/user/Desktop/variable.Rdata")

binary_interaction <- sum_interaction
save(bone_interaction, brain_interaction, kidney_interaction, liver_interaction, lung_interaction,
     binary_interaction, count_interaction, file = "C:/Users/user/Desktop/interaction.Rdata")
















