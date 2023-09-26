file_name <- file.choose() # Choose the original data set
data <- readRDS(file_name) # and name "data"

head(data)
tail(data)

attach(data)
uniq_zip <- unique(zip) # Find unique zip codes
num_zip_rep <- rep(NA, length(uniq_zip))
data_prcd <- matrix(NA, nrow = length(uniq_zip), ncol=25) # Make a matrix (zip x variable) to store the processed data

# aggregate the original data to zip code level
for(zz in 1:length(uniq_zip)){
  temp_subset <- data[data$zip==unique(zip)[zz],] # Find observations with the same zip code
  data_prcd[zz,] <- c(
    # Outcome: death rate = (dead)/(total person year)
    death_rate = sum(temp_subset$dead)/sum(temp_subset$time_count), 
    
    # Covariates 
    # NOTE: These are sample means of each covariate which would not represent true values.
    #       It is more appropriate to calculate weighted averages, taking into account the number of individuals in each original unit.
    pct_female = mean(temp_subset$sex-1), 
    pct_dual = mean(temp_subset$dual), # clarify it's okay to disregard entry age, race 
    t(apply(temp_subset[12:25],2,mean)),
    
    # Exposures (Multivariate Continuous)
    # NOTE: We also need to figure out how to define exposures for each zip code.
    t(apply(temp_subset[11],2,mean)), # pm2.5 ensemble
    exposure2 = rnorm(1), # don't want to do this anymore, want to use value of components 
    # t(apply(temp_subset[25::],2,mean)),
    exposure3 = rnorm(1),
    exposure4 = rnorm(1),
    exposure5 = rnorm(1),
    exposure6 = rnorm(1),
    exposure7 = rnorm(1),
    exposure8 = rnorm(1))
  cat(zz, "\r")
}

colnames(data_prcd) <- c("death_rate", "pct_female", "pct_dual",
                         dimnames(data[12:25])[[2]],
                         dimnames(data[11])[[2]],
                         paste0("exposure", 2:8))

