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










