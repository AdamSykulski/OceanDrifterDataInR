# ACTIVITY 1: COMPARING DRIFTER AND SATELLITE SEA SURFACE TEMPERATURE (SST) DATA

# The first goal of this activity is to compare two SST estimate products:
# - The Global Drifter Program hourly SST product (Elipot et al. 2022)
# - the Multi-scale Ultra-high Resolution (MUR) SST analysis (Chin et al. 2017)
# Each product also reports associated measures of uncertainty (given as standard deviations of the estimates)
# The second goal of this activity is to interpret these uncertainties and use them in the comparison.

# Learning Outcomes
# After this computer lab session you should be able to
# - perform statistical comparisons of two datasets in R measuring the same quantity
# - interpret uncertainty estimates and use them to establish if datasets are significantly different from each other

# Load some libraries (uncomment the installation lines if you have not installed these packages before)
# install.packages("ncd4")
library(ncdf4)
# install.packages("ggplot2")
library(ggplot2)

# Load data, which correspond to 74921 SST measurements in the North Atlantic
sst <- ncvar_get(nc_open("gdp_9am_atlantic.nc"), "sst") # SST from drifters
err_sst <- ncvar_get(nc_open("gdp_9am_atlantic.nc"), "err_sst") # corresponding uncertainties from drifters
mur_sst <- ncvar_get(nc_open("gdp_9am_atlantic.nc"), "mur_sst") # corresponding SST from satellites
mur_sst_err <- ncvar_get(nc_open("gdp_9am_atlantic.nc"), "mur_sst_err") # corresponding uncertainties from satellites

# Preprocessing: Remove NaNs (from all) and outliers (>2K difference)
a <- is.na(err_sst)
sst <- sst[!a,drop = FALSE]
err_sst <- err_sst[!a,drop = FALSE]
mur_sst_err <- mur_sst_err[!a,drop = FALSE]
mur_sst <- mur_sst[!a,drop = FALSE]
A <- mur_sst - sst # the difference
CT <- 2 # cutoff difference in K
sst <- sst[abs(A) < CT]
err_sst <- err_sst[abs(A) < CT]
mur_sst_err <- mur_sst_err[abs(A) < CT]
mur_sst <- mur_sst[abs(A) < CT]
proportion_data_left <- sum(abs(A) < CT) / length(A) # with 2K cuttoff, 99.3% data retained

# Q1: Visualise the differences between the drifter and satellite SST data.
# As examples, plot comparisons such as scatter plots (add a y=x line), 2-d histograms,
# and a histogram of the differences overlaying a best-fit normal distribution.
# What initial conclusions do you have?



# Q2: What is the average and standard deviation of the differences between the observations?
# Perform a hypothesis test to establish whether the average difference is significantly different from zero



# Q3: Now compute the average reported uncertainties of each product, which product seems to have higher uncertainty?
# Also plot histograms of their distributions (suggest to plot drifter uncertainties on a dB scale)
# Finally, plot a scatter plot or 2D histogram of the uncertainties against each other, is an independence assumption reasonable?



# Q4: Assuming independence, what is the associated uncertainty of the differences between observations? Consult your lecturer if unsure!
# are these uncertainties broadly consistent with the observed standard deviations of the differences?
# what might you do next to establish if the reported uncertainties and observed standard deviations agree?


