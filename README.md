## Data Processing
We used `data_processing.R` to process the data so that the data structure is similar to that of
- [Ren, B., Wu, X., Braun, D., Pillai, N., & Dominici, F. (2021). Bayesian modeling for exposure response curve via Gaussian processes: causal effects of exposure to air pollution on health outcomes. arXiv preprint arXiv:2105.03454.](https://arxiv.org/pdf/2105.03454.pdf)

We aggregate the data to the zip code level by taking the sample mean of all observations sharing the same zip code **for testing the code at this time**.

![image](https://github.com/hshin111/Medicare_Analysis_Code_1/assets/56696926/d6823cc9-db5a-4a0c-9bf6-743edc9d4a1f)
For example, there are 19 units in zip code 40410 and we define the outcome (death rate) by a weighted average of "dead"s where weights are given by "time_count" whereas covariates and exposures are calculated by taking sample mean without weights. Then, the information about zip code 40410 is transformed as

![image](https://github.com/hshin111/Medicare_Analysis_Code_1/assets/56696926/bcf9c344-14a5-4c1f-b091-33f230f46b1e)

We also added randomly generated 7 continuous exposures to the data so that we have 1 outcome variable, 16 covariates, and 8 exposures all continuous. After the processing, we have 25,510 units and we expect this would not exceed 50k as the number of zip codes in the US is less than that.

## SepBART Function
`00-functions-block` contains the main function. To run the code, `ModBart` and `SoftBart` (MFM version, https://github.com/theodds/SoftBART/tree/MFM) of Antonio Linero, and `Rcpp` are required. An example usage can be found in `sim.R` file.


# Processed Data Example with randomized values 

**Data Table Format Explanation**

In this readme, we provide an explanation regarding the format of the data table. Please note that the table format described here is in the long version, as opposed to the wide version that was used during the data processing phase.

**Wide vs. Long Format**

- **Wide Format**: During data processing, we utilized the wide format of the table. In the wide format, each row typically represents one observation, and each column represents a distinct variable. This format is often convenient for initial data manipulation and calculations.

- **Long Format**: However, for the sake of clarity and organization, we have transformed the data table into the long format for presentation and documentation purposes. In the long format, data is organized differently: it is structured to have a column for variable names, a column for corresponding values, and additional columns for identifying information, such as year and location. This format is especially useful when displaying or sharing data for analysis or reporting.

**Year-Specific Datasets**

It's important to note that the data presented here pertains to a single year. Our dataset is organized such that each year is treated as a separate dataset. This separation by year allows for easy access and analysis of data for specific time periods. As such, the long format data table provided here simplifies the process of accessing and interpreting information for individual years.


| Variable                                    | zip_1  | zip_2  | zip_3  |
|---------------------------------------------|--------|--------|--------|
| pm25_ensemble                               | 12.5   | 15.2   | 10.8   |
| mean_bmi                                    | 26.7   | 28.1   | 25.4   |
| smoke_rate                                  | 18.3   | 14.7   | 21.5   |
| hispanic                                    | 8.2    | 12.6   | 5.7    |
| pct_blk                                     | 15.8   | 22.3   | 10.1   |
| medhouseholdincome                          | 75000  | 60000  | 80000  |
| medianhousevalue                            | 450000 | 350000 | 550000 |
| poverty                                     | 12.5   | 20.3   | 8.7    |
| education                                   | 14.5   | 12.8   | 16.3   |
| popdensity                                  | 1200   | 800    | 1500   |
| pct_owner_occ                               | 65.2   | 45.7   | 72.1   |
| summer_tmmx                                | 89.4   | 92.0   | 86.8   |
| winter_tmmx                                | 65.1   | 58.7   | 70.3   |
| summer_rmax                                | 4.5    | 5.2    | 4.1    |
| winter_rmax                                | 3.2    | 2.8    | 3.6    |
| pm25                                       | 18.7   | 22.3   | 15.9   |
| ozone                                      | 0.045  | 0.055  | 0.038  |
| br                                         | 0.002  | 0.001  | 0.003  |
| ca                                         | 0.005  | 0.004  | 0.006  |
| cu                                         | 0.0005 | 0.0004 | 0.0006 |
| ec                                         | 0.0012 | 0.0015 | 0.0010 |
| fe                                         | 0.0023 | 0.0021 | 0.0025 |
| k                                          | 0.0010 | 0.0008 | 0.0012 |
| nh4                                        | 0.0035 | 0.0032 | 0.0037 |
| ni                                         | 0.0003 | 0.0002 | 0.0004 |
| no3                                        | 0.006  | 0.005  | 0.007  |
| oc                                         | 0.0055 | 0.0062 | 0.0048 |
| pb                                         | 0.0001 | 0.00009| 0.00012|
| si                                         | 0.003  | 0.0027 | 0.0033 |
| so4                                        | 0.008  | 0.0075 | 0.0083 |
| v                                          | 0.0004 | 0.0003 | 0.0005 |
| z                                          | 0.0032 | 0.0038 | 0.0029 |
| randall_martin_PM25                        | 14.2   | 16.8   | 13.4   |
| randall_martin_BC                          | 0.0025 | 0.0020 | 0.0028 |
| randall_martin_NH4                         | 0.0040 | 0.0036 | 0.0043 |
| randall_martin_NIT                         | 0.007  | 0.0063 | 0.0076 |
| randall_martin_OM                          | 0.0052 | 0.0057 | 0.0048 |
| randall_martin_SO4                         | 0.0092 | 0.0085 | 0.0097 |
| randall_martin_SOIL                        | 0.003  | 0.0028 | 0.0032 |
| randall_martin_SS                          | 0.0015 | 0.0013 | 0.0016 |
| total_count                                | 850    | 1200   | 750    |
| female_percentage                          | 52.4   | 45.8   | 53.2   |
| dual_percentage                            | 18.7   | 22.5   | 17.9   |
| mean_age                                   | 72.3   | 68.7   | 73.8   |
| death_rate                                 | 10.5   | 15.2   | 9.8    |
| percentage_race_labelAsian                 | 5.2    | 8.1    | 4.3    |
| percentage_race_labelBlack                 | 12.8   | 17.3   | 11.1   |
| percentage_race_labelOther                 | 3.7    | 2.9    | 4.5    |
| percentage_race_labelWhite                 | 77.9   | 69.7   | 79.2   |
| percentage_race_labelHispanic              | 8.3    | 11.4   | 7.4    |
| percentage_race_labelNorth_American_Native | 0.1    | 0.4    | 0.2    |


# Processed Data Variable Descriptions

## Demographic and Socioeconomic Variables

- `zip`: Zip code.
- `year`: Year of data collection.
- `pm25_ensemble`: Ensemble PM2.5 data.
- `mean_bmi`: Mean Body Mass Index (BMI) of the population.
- `smoke_rate`: Rate of smoking in the population.
- `hispanic`: Percentage of the population identifying as Hispanic.
- `pct_blk`: Percentage of the population identifying as Black.
- `medhouseholdincome`: Median household income.
- `medianhousevalue`: Median value of houses in the area.
- `poverty`: Poverty rate in the population.
- `education`: Education level of the population.
- `popdensity`: Population density.
- `pct_owner_occ`: Percentage of owner-occupied households.
- `mean_age`: Mean age of the population.
- `female_percentage`: Percentage of the population identifying as female.
- `dual_percentage`: Percentage of individuals with dual proxy to poverty.

## Environmental Variables

- `summer_tmmx`: Maximum summer temperature.
- `winter_tmmx`: Maximum winter temperature.
- `summer_rmax`: Maximum summer precipitation.
- `winter_rmax`: Maximum winter precipitation.
- `pm25`: Particulate Matter (PM2.5) concentration.
- `ozone`: Ozone concentration.
- `br`, `ca`, `cu`, `ec`, `fe`, `k`, `nh4`, `ni`, `no3`, `oc`, `pb`, `si`, `so4`, `v`, `z`: Concentrations of different chemical elements with Joel Schwartz data.

## Randall Martin Variables

- `randall_martin_PM25`: PM2.5 data from the Randall Martin data.
- `randall_martin_BC`: Black Carbon (BC) data from the Randall Martin data.
- `randall_martin_NH4`: Ammonium (NH4) data from the Randall Martin data.
- `randall_martin_NIT`: Nitrate (NIT) data from the Randall Martin data.
- `randall_martin_OM`: Organic Matter (OM) data from the Randall Martin data.
- `randall_martin_SO4`: Sulfate (SO4) data from the Randall Martin data.
- `randall_martin_SOIL`: Soil data from the Randall Martin data.
- `randall_martin_SS`: Suspended Soil (SS) data from the Randall Martin data.

## Health and Mortality Variables

- `total_count`: Total count of individuals.
- `death_rate`: Death rate in the population.
- `percentage_race_labelAsian`: Percentage of the population identifying as Asian.
- `percentage_race_labelBlack`: Percentage of the population identifying as Black.
- `percentage_race_labelOther`: Percentage of the population identifying as Other race.
- `percentage_race_labelWhite`: Percentage of the population identifying as White.
- `percentage_race_labelHispanic`: Percentage of the population identifying as Hispanic.
- `percentage_race_labelNorth_American_Native`: Percentage of the population identifying as North American Native.



# Data Processing Steps


## Step 1: Load and Preprocess Medicare Data
```R
# Load the Medicare data from an RDS file
file_path <- "../data/input/aggregated_2000-2016_medicare_mortality_pm25_zip/aggregate_data.RDS"
data <- readRDS(file_path)

# Additional preprocessing steps if required.
```

## Step 2: Load and Preprocess PM2.5 Data
```R
# Initialize an empty list to store data frames for each year
data_list <- list()

# Loop through years from 2000 to 2016
for (year in 2000:2016) {
  # File path to the RDS file for the current year
  file_path <- paste0("../data/input/PM25_v2/annual/", year, ".rds")
  
  # Read the data from the RDS file and preprocess
  # Additional preprocessing steps if required.
  
  # Append the data frame to the list
  data_list[[year]] <- annual_pm25
}

# Combine all data frames into a single data frame named 'joel_schwartz_pm25'
joel_schwartz_pm25 <- do.call(rbind, data_list)
```

## Step 3: Load and Preprocess PM2.5 Components Data
```R
# Similar steps as in Step 2 but for PM2.5 components data
# Combine all data frames into a single data frame named 'joel_schwartz_pm25_components'
```

## Step 4: Load and Preprocess Annual O3 Data
```R
# Similar steps as in Step 2 but for O3 data
# Combine all data frames into a single data frame named 'joel_schwartz_annual_O3'
```

## Step 5: Load and Preprocess Randall Martin PM2.5 Components Data
```R
# Similar steps as in Step 2 but for Randall Martin PM2.5 components data
# Combine all data frames into a single data frame
```

## Step 6: Merging Data by Year
```R
# Load and preprocess data for each year
# Merge the data frames using the zip code as the key
# Additional data cleaning and preprocessing steps
# Save the merged data for each year as an RDS file
```

## Step 7: Final Data Merge
```R
# Merge all the yearly data into a single dataset for the entire period
# Additional data cleaning and preprocessing steps
# Save the final merged data as an RDS file
```

# SepBART Processing Decision

1. **Y (Outcome Variable):**
   - **Y <- data_prcd[,49] # death rate (percentage)**
   - Y represents the outcome variable, which is the "death rate" expressed as a percentage. This variable likely contains data related to mortality rates.

2. **X (Covariates):**
   - **X <- data_prcd[, c(46:48, 50:55, 4:5, 8:13)]**
   - **X_numeric <- matrix(as.numeric(X), nrow = nrow(X))**
   - X represents a set of covariates or independent variables used in statistical analysis. These covariates are numeric in nature and include various factors, such as "mean_bmi," "smoke_rate," "medhouseholdincome," "medianhousevalue," "poverty," "education," "popdensity," "pct_owner_occ," "female_percentage," "dual_percentage," "mean_age," and demographic information like "percentage_race_labelAsian," "percentage_race_labelBlack," "percentage_race_labelOther," "percentage_race_labelWhite," "percentage_race_labelHispanic," and "percentage_race_labelNorth_American_Native."

3. **W (Exposure Variables):**
   - **W <- data_prcd[, c(14:34, 37:44)]**
   - W represents exposure variables. These variables are often used to examine how different factors or exposures influence the outcome variable. In this case, the exposure variables include "summer_tmmx," "winter_tmmx," "summer_rmax," "winter_rmax," "pm25," "ozone," "br," "ca," "cu," "ec," "fe," "k," "nh4," "ni," "no3," "oc," "pb," "si," "so4," "v," "z," "randall_martin_PM25," "randall_martin_BC," "randall_martin_NH4," "randall_martin_NIT," "randall_martin_OM," "randall_martin_SO4," "randall_martin_SOIL," and "randall_martin_SS." These variables may represent environmental, chemical, or other exposures of interest.

In summary, X represents the covariates or independent variables, Y represents the outcome variable, and W represents exposure variables. These variables are often used in statistical modeling and analysis to understand relationships between different factors and the outcome of interest, which, in this case, is the death rate expressed as a percentage.
