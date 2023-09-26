library("ModBart")
library("SoftBart") # MFM version
library("mvtnorm")
library("Rcpp")

## Call the SepBART function
source('../00-functions-block.R')

## Load the data
data <- readRDS("../data/intermediate/SepBART_data/preprocessed_data/all_merged_2000.rds")
cat("fitting data for year 2000")
# Generate 1000 random row indices
#random_indices <- sample(nrow(data), size = 100, replace = FALSE)

# Subset the dataframe using the random indices
#random_sample <- data[random_indices, ]

#data_prcd <- as.matrix(random_sample)

data_prcd <- as.matrix(data)

Y <- data_prcd[,49] # death rate (percentage)

X <- data_prcd[, c(46:48,50:55, 4:5, 8:17)]
# Convert character matrix to a numeric matrix element-wise
X_numeric <- matrix(as.numeric(X), nrow = nrow(X))

W <- data_prcd[,c(19:25,27,29:31,33:34, 37:44)]
W_numeric <- matrix(as.numeric(W), nrow = nrow(W))

W0 <- apply(W_numeric, 2, quantile, 0.4, na.rm = TRUE)
W1 <- apply(W_numeric, 2, quantile, 0.6, na.rm = TRUE)

Y_numeric <- as.numeric(Y)
X_numeric <- matrix(as.numeric(X), nrow = nrow(X))
W1_numeric <- as.numeric(W1)
W0_numeric <- as.numeric(W0)

## Fit the model.
fit <- SepBART_blk(Y = Y_numeric, X = X_numeric, W = W_numeric,
                   Xtest = NULL, W1 = W1_numeric, W0 = W0_numeric,
                   nMCMC=1000,BurnIn_portion=0.2,stepsize=5)


# Define the file path to save the 'fit' object
output_directory <- "../data/output/"
dir.create(output_directory, recursive = TRUE)  # Create the directory if it doesn't exist
output_file_path <- paste0(output_directory, "fit_result_2000.rds")

# Save the 'fit' object as an RDS file
saveRDS(fit, output_file_path)

# Print a message to confirm successful saving
cat("The 'fit' object has been saved to:", output_file_path, "\n")
