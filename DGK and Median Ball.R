# Load necessary libraries
library(MASS) # for the ginv() function
library(matrixStats)  # Untuk fungsi rowMedians
library(readxl)
library(ggplot2)
library(GGally)
library(tidyr)
library(dplyr)
library(psych)
library(MVN)
library(reshape2)
library(rrcov)

# Import Data
df <- read_excel("PDRB atas Dasar Harga Baku Menurut Pengeluaran (Juta Rupiah) per Provinsi Tahun 2023.xlsx",
                 sheet = "Sheet2")

options(digits = 2)

data=df[, 2:8]

data

#Analisis Deskriptif
options(scipen = 999)
summarydata<-summary(df[,-1])
summarydata<-as.data.frame(summarydata)
summarydata
write.csv(summarydata, "summarydata.csv", row.names = FALSE)
options(scipen=0, digits=8)

# Membuat vektor nama baru untuk variabel
new_names <- c("PK-RT", "PK-LNPRT", "PK-P", "PMTB", "PI", "NetEkspor", "PDRB")

# Fungsi panel untuk menambahkan garis linear
panel.smooth <- function(x, y) {
  points(x, y)
  abline(lm(y ~ x), col = "red", lwd = 2)  # Menambahkan garis regresi linear
}

# Membuat scatter plot matrix dengan garis linear
pairs(data, labels = new_names, main = "Scatter Plot Matrix Variabel PDRB Pengeluaran", panel = panel.smooth)

# Membuat boxplot untuk setiap variabel
boxplot(data, main="Boxplot Setiap Variabel", 
        las=1, col="lightblue",names=new_names)

# Function to compute Mahalanobis distance
mahalanobis_distance <- function(X, mean_vec, cov_mat) {
  # Ensure inputs are numeric
  X <- as.matrix(X)
  mean_vec <- as.numeric(mean_vec)
  cov_mat <- as.matrix(cov_mat)
  
  diff <- sweep(X, 2, mean_vec)
  dist <- sqrt(rowSums((diff %*% ginv(cov_mat)) * diff))
  return(dist)
}

#Check Outlier
m <- as.matrix(data)
m

# Initial estimation
cm <- colMeans(m)
cm
cv <- cov(m)
cv

# Mahalanobis distance
mahd<-mahalanobis_distance(m,cm,cv)
mahd

# Menambahkan jarak Mahalanobis sebagai kolom baru
data_with_md <- cbind(data, Mahalanobis_Distance = mahd)
data_with_md

# Menentukan threshold menggunakan distribusi chi-square
alpha <- 0.95
threshold <- qchisq(1-alpha, df = ncol(data))
threshold

# Mengidentifikasi outlier
outliers <- mahd > threshold

# Menambahkan keterangan "outlier" atau "bukan outlier"
outlier_labels <- ifelse(outliers, "Outlier", "Bukan outlier")

# Menambahkan informasi outlier ke dalam data
data_with_md <- cbind(data_with_md, Outlier = outliers)
data_with_md

write.csv(data_with_md, "data_with_md1.csv", row.names = FALSE)

# Membuat data frame untuk plotting
provinsi<-df[,1]
data_with_md2 <- data.frame(
  Provinsi = provinsi,
  Observation = 1:nrow(data),
  Mahalanobis_Distance = mahd,
  Outlier = factor(outliers, levels = c(FALSE, TRUE), labels = c("FALSE", "TRUE"))
)

# Reset dan buka perangkat grafik baru
dev.off() 
windows()

# Plot sebaran jarak Mahalanobis
ggplot(data_with_md2, aes(x = Observation, y = Mahalanobis_Distance, color = Outlier)) +
  geom_point(size = 5) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
  geom_text(aes(label = Provinsi), hjust = -0.1, vjust = 0, size = 3) + 
  labs(
    title = "Plot of Mahalanobis Distances",
    x = "Provinsi",
    y = "Mahalanobis Distance"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

ggplot(data_with_md2, aes(x = Observation, y = Mahalanobis_Distance, color = Outlier)) +
  geom_point(size = 5) +
  geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
  geom_text(
    aes(label = ifelse(Outlier == TRUE, as.character(Provinsi), "")), 
    hjust = -0.1, 
    vjust = 0, 
    size = 5
  ) +
  labs(
    title = "Plot of Mahalanobis Distances",
    x = "Provinsi",
    y = "Mahalanobis Distance"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

# DGK Algorithm
dgk_algorithm <- function(X, tol = 1e-6, max_iter = 100) {
  X <- as.matrix(X)  
  n <- nrow(X)
  p <- ncol(X)
  
  # Initial estimation
  T_start <- colMeans(X)
  C_start <- cov(X)
  
  T_k <- T_start
  C_k <- C_start
  
  mahalanobis_distances <- list()  
  covariance_matrices <- list()    
  converged <- FALSE
  
  for (k in 1:max_iter) {
    # Compute Mahalanobis distance
    D_k <- mahalanobis_distance(X, T_k, C_k)
    mahalanobis_distances[[k]] <- D_k  
    covariance_matrices[[k]] <- C_k    
    
    # Compute median of Mahalanobis distance
    Med_k <- median(D_k)
    
    # Select new data set
    X_new <- X[D_k <= Med_k, ]
    
    # Update mean and covariance
    T_new <- colMeans(X_new)
    C_new <- cov(X_new)
    
    # Check for convergence
    if (all(abs(diag(C_new) - diag(C_k)) < tol)) {
      converged <- TRUE
      covariance_matrices[[k + 1]] <- C_new  
      break
    }
    
    # Update T_k and C_k
    T_k <- T_new
    C_k <- C_new
    
    # Check if diagonal of C_new equals diagonal of C_start
    if (all(abs(diag(C_new) - diag(C_start)) < tol)) {
      converged <- TRUE
      covariance_matrices[[k + 1]] <- C_new  # Store the final covariance matrix
      break
    }
  }
  
  if (!converged) {
    warning("Algorithm did not converge within the maximum number of iterations")
  }
  
  return(list(mean = T_k, covariance = C_k, iterations = k, 
              mahalanobis_distances = mahalanobis_distances, 
              covariance_matrices = covariance_matrices))
}

# Apply DGK Algorithm
result <- dgk_algorithm(data)

# Print means and covariances at the final iteration
cat("Final Mean Vector:\n")
print(result$mean)
cat("\nFinal Covariance Matrix:\n")
print(result$covariance)
cat("\nFinal Mahalanobis Distance:\n")
print(result$mahalanobis_distances)

# Print covariance matrices for each iteration
for (i in 1:length(result$covariance_matrices)) {
  cat("\nCovariance Matrix at Iteration", i, ":\n")
  print(result$covariance_matrices[[i]])
}

# Print covariance matrices for each iteration
for (i in 1:length(result$mean)) {
  cat("\nmean at Iteration", i, ":\n")
  print(result$mean[[i]])
}

# Create a data frame from the original matrix
X_with_distances <- as.data.frame(data)

# Iterations
num_iterations <- length(result$mahalanobis_distances)

# Add MD to Dataset
for (i in 1:num_iterations) {
  distance_column_name <- paste("Mahalanobis_Distance_Iter", i, sep = "_")
  X_with_distances[[distance_column_name]] <- result$mahalanobis_distances[[i]]
}

# Save DGK with MD
print(head(X_with_distances))
write.csv(X_with_distances, "X_with_distances.csv", row.names = FALSE)


# MB Algorithm
mb_algorithm <- function(X, tol = 1e-6, max_iter = 100) {
  X <- as.matrix(X) 
  n <- nrow(X)
  p <- ncol(X)
  
  # Matrix Identitas dan Median
  T_start <- apply(X, 2, median)
  C_start <- diag(p)
  
  T_k <- T_start
  C_k <- C_start
  
  mahalanobis_distances <- list()  
  covariance_matrices <- list()    
  converged <- FALSE
  
  for (k in 1:max_iter) {
    # Compute Mahalanobis distance
    D_k <- mahalanobis_distance(X, T_k, C_k)
    mahalanobis_distances[[k]] <- D_k  
    covariance_matrices[[k]] <- C_k    
    
    # Compute median of Mahalanobis distance
    Med_k <- median(D_k)
    
    # Select new data set
    X_new <- X[D_k <= Med_k, ]
    
    # Update mean and covariance
    T_new <- colMeans(X_new)
    C_new <- cov(X_new)
    
    # Check for convergence
    if (all(abs(diag(C_new) - diag(C_k)) < tol)) {
      converged <- TRUE
      covariance_matrices[[k + 1]] <- C_new  
      break
    }
    
    # Update T_k and C_k
    T_k <- T_new
    C_k <- C_new
    
    # Check if diagonal of C_new equals diagonal of C_start
    if (all(abs(diag(C_new) - diag(C_start)) < tol)) {
      converged <- TRUE
      covariance_matrices[[k + 1]] <- C_new 
      break
    }
  }
  
  if (!converged) {
    warning("Algorithm did not converge within the maximum number of iterations")
  }
  
  return(list(mean = T_k, covariance = C_k, iterations = k, 
              mahalanobis_distances = mahalanobis_distances, 
              covariance_matrices = covariance_matrices))
}

# Apply MB Algorithm on data
resultMB <- mb_algorithm(data)

# Print means and covariances at the final iteration
cat("Final Mean Vector:\n")
print(resultMB$mean)
cat("\nFinal Covariance Matrix:\n")
print(resultMB$covariance)
cat("\nFinal Mahalanobis Distance:\n")
print(resultMB$mahalanobis_distances)

# Print covariance matrices for each iteration
for (i in 1:length(resultMB$covariance_matrices)) {
  cat("\nCovariance Matrix at Iteration", i, ":\n")
  print(resultMB$covariance_matrices[[i]])
}

# Calculate the determinant of the final covariance matrix
final_cov_matrix <- resultMB$covariance
det_final_cov_matrix <- det(final_cov_matrix)
cat("\nDeterminant of Final Covariance Matrix:\n")
print(sqrt(det_final_cov_matrix))

# Create a dataframe from the original matrix
XMB_with_distances <- as.data.frame(data)

# Ensure there are enough iterations
num_iterations <- length(resultMB$mahalanobis_distances)

# Add Mahalanobis distance for each iteration as new columns
for (i in 1:num_iterations) {
  distance_column_name <- paste("Mahalanobis_Distance_Iter", i, sep = "_")
  XMB_with_distances[[distance_column_name]] <- resultMB$mahalanobis_distances[[i]]
}

# Print the resulting dataset with Mahalanobis distance for each iteration
print(head(XMB_with_distances))
write.csv(XMB_with_distances, "XMB_with_distances.csv", row.names = FALSE)

#RFCH Algorithm
# Define a function for RFCH algorithm
rfch_algorithm <- function(X, DGK_result, MB_result, tol = 1e-6, max_iter = 100) {
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  # Calculate determinant for DGK and MB
  det_DGK <- det(DGK_result$covariance)
  det_MB <- det(MB_result$covariance)
  
  # Initial estimator FCH
  if (det_DGK < det_MB) {
    T_0_FCH <- DGK_result$mean
    C_0_FCH <- (median(DGK_result$mahalanobis_distances[[DGK_result$iterations]]) / qchisq(0.5, p)) * DGK_result$covariance
  } else {
    T_0_FCH <- MB_result$mean
    C_0_FCH <- (median(MB_result$mahalanobis_distances[[MB_result$iterations]]) / qchisq(0.5, p)) * MB_result$covariance
  }
  
  T_k_FCH <- T_0_FCH
  C_k_FCH <- C_0_FCH
  C_k_RFCH <- C_0_FCH
  mahalanobis_distances <- list()
  covariance_matrices <- list()
  
  converged <- FALSE
  k <- 0
  
  while (!converged && k < max_iter) {
    k <- k + 1
    
    # Compute Mahalanobis distance
    D_k_FCH <- mahalanobis_distance(X, T_k_FCH, C_k_FCH)
    mahalanobis_distances[[k]] <- D_k_FCH
    covariance_matrices[[k]] <- C_k_FCH
    
    # Selection of data based on chi-squared threshold
    threshold <- qchisq(0.975, p)
    X_new <- X[D_k_FCH <= threshold, , drop = FALSE]
    
    # Update mean and covariance
    T_k_FCH <- colMeans(X_new)
    C_k_FCH <- cov(X_new)
    
    # Reweighting step
    C_k_RFCH_new <- (median(D_k_FCH) / qchisq(0.5, p)) * C_k_FCH
    
    # Check for convergence
    if (all(abs(diag(C_k_RFCH_new) - diag(C_k_RFCH)) < tol)) {
      converged <- TRUE
    } else {
      C_k_RFCH <- C_k_RFCH_new
    }
  }
  
  if (!converged) {
    warning("RFCH algorithm did not converge within the maximum number of iterations.")
  }
  
  return(list(mean = T_k_FCH, covariance = C_k_RFCH, iterations = k, 
              mahalanobis_distances = mahalanobis_distances, 
              covariance_matrices = covariance_matrices,
              converged = converged,
              X_final = X_new))
}

# Apply DGK Algorithm
dgk_result <- dgk_algorithm(data)
str(dgk_result)

# Apply MB Algorithm
mb_result <- mb_algorithm(data)
str(mb_result)

# Apply RFCH Algorithm
rfch_result <- rfch_algorithm(data, dgk_result, mb_result)

# Print means and covariances at the final iteration
cat("Final Mean Vector for RFCH:\n")
print(rfch_result$mean)
cat("\nFinal Covariance Matrix for RFCH:\n")
print(rfch_result$covariance_matrices[[6]])

# Print covariance matrices for each iteration
for (i in 1:length(rfch_result$covariance_matrices)) {
  cat("\nCovariance Matrix at Iteration", i, ":\n")
  print(rfch_result$covariance_matrices[[i]])
}

# Create a dataframe from the original matrix
X_with_distances_rfch <- as.data.frame(data)

# Ensure there are enough iterations
num_iterations_rfch <- length(rfch_result$mahalanobis_distances)

# Add Mahalanobis distance for each iteration as new columns
for (i in 1:num_iterations_rfch) {
  distance_column_name <- paste("Mahalanobis_Distance_RFCH_Iter", i, sep = "_")
  X_with_distances_rfch[[distance_column_name]] <- rfch_result$mahalanobis_distances[[i]]
}

# Print the resulting dataset with Mahalanobis distance for each iteration
print(head(X_with_distances_rfch))
write.csv(X_with_distances_rfch, "X_with_distancesRFCH.csv", row.names = FALSE)

# save the dataframe to a CSV file
newww<-as.data.frame(rfch_result$covariance_matrices[[6]])
write.csv(newww, "X_with_distances_rfch.csv", row.names = FALSE)

# Membuat scatter plot matrix dengan garis linear
pairs(finalRFCH, labels = new_names, main = "Scatter Plot Matrix Variabel PDRB Pengeluaran", panel = panel.smooth)


# Define the correlation matrix
corrr_matrix <- matrix(c(
  1, 0.969, 0.915, 0.954, 0.828, 0.292, 0.975,
  0.969, 1, 0.937, 0.912, 0.812, 0.22, 0.934,
  0.915, 0.937, 1, 0.906, 0.839, 0.283, 0.912,
  0.954, 0.912, 0.906, 1, 0.779, 0.504, 0.991,
  0.828, 0.812, 0.839, 0.779, 1, 0.233, 0.81,
  0.292, 0.22, 0.283, 0.504, 0.233, 1, 0.489,
  0.975, 0.934, 0.912, 0.991, 0.81, 0.489, 1
), nrow = 7, byrow = TRUE)

# Name the rows and columns
variable_names=c("X1", "X2",
                 "X3", "X4", 
                 "X5", "X6", 
                 "X7")
rownames(corrr_matrix) <- colnames(corrr_matrix) <- variable_names

# Reshape the data using melt from reshape2
cordat <- melt(corrr_matrix)

# Create the heatmap using ggplot2
heatmap <- ggplot(data = cordat, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 4, color = "black")+
  scale_fill_gradient2(low = "#00008B", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0,1), space = "Lab", 
                       name="Correlation") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(vjust = 1, size = 12, hjust = 1),
    axis.text.y = element_text(vjust = 1, size = 12, hjust = 1),
    axis.title.x = element_blank(), 
    axis.title.y = element_blank()
  )  +
  coord_fixed()

# Print the heatmap
print(heatmap)

# Check Distribution
finalRFCH<-rfch_result$X_final
head(finalRFCH)
mvn(finalRFCH, mvnTest = "mardia")
head(df)
mvn(data, mvnTest = "mardia")

# Menghitung korelasi Pearson
cor_matrix <- cor(data, method = "pearson")
cor_matrixx<-as.data.frame(cor_matrix)
write.csv(cor_matrix, "corrrr.csv", row.names = FALSE)

# Mencetak matriks korelasi
print(cor_matrix)

# Untuk visualisasi yang lebih baik, kita bisa menggunakan library corrplot
library(corrplot)

# Membuat plot korelasi
corrplot(cor_matrix, method = "color", type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
