# ACTIVITY 1: COMPARING DRIFTER AND SATELLITE SEA SURFACE TEMPERATURE (SST) DATA

# The first goal of this activity is to compare two SST estimate products:
# - The Global Drifter Program hourly SST product (Elipot et al. 2022)
# - the Multi-scale Ultra-high Resolution (MUR) SST analysis (Chin et al. 2017)
# Each product also reports associated measures of uncertainty (given as standard errors of the estimates)
# The second goal of this activity is to interpret these uncertainties and use them in the comparison.

# Learning Outcomes
# After this computer lab session you should be able to
# - perform statistical comparisons of two datasets in R measuring the same quantity
# - interpret uncertainty estimates and use them to establish if datasets are significantly different from each other

# Load some libraries (uncomment the installation lines if you have not installed these packages before)
# install.packages("ncdf4")
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
plot(sst, mur_sst, xlab='drifter', ylab='MUR') # scatter plot
abline(a=0, b=1, col='red', lwd=2) # add y=x line
data <- data.frame( x=sst, y=mur_sst ) # data frame for 2-d histogram
ggplot(data, aes(x=sst, y=mur_sst) ) +
  geom_bin2d(bins = 200) +
  scale_fill_continuous(type = "viridis") +
  theme_bw() # 2-d histogram in ggplot
diff = sst-mur_sst # differences between measurements
hist(diff, breaks=100, probability=TRUE, main='', xlab='MUR-drifter') # histogram
curve(dnorm(x, mean(diff), sd(diff)), add=TRUE, col='blue', lwd=2) # add the normal distribution

# Q2: What is the average and standard deviation of the differences between the observations?
# Perform a hypothesis test to establish whether the average difference is significantly different from zero
obs_difference <- mean(diff) # average observed differences
obs_standard_error <- sd(diff) # average observed standard deviation
sem <- obs_standard_error / sqrt(length(diff)) # standard error of the mean to perform a z-test (could do a t-test too, same answer)
pval <- 2*pnorm(-abs(obs_difference) / sem) # p-value < 0.05 so reject, there is a (small) statistically significant difference
# the satellite estimates seem to be reporting slightly lower SST than drifters on average

# Q3: Now compute the average reported standard error of each product, which product seems to have the larger errors?
# Also plot histograms of their distributions (suggest to plot drifter errors on a dB scale)
# Finally, plot a scatter plot or 2D histogram of the errors against each other, is an independence assumption reasonable?
mur_average_standard_error <- mean(mur_sst_err) # slightly lower than observed standard error
drifter_average_standard_error <- mean(err_sst) # very low!
hist(mur_sst_err, main='', xlab='MUR errors') # histograms
hist(10 * log10(err_sst), main='', xlab='drifter errors (dB)')
plot(10*log10(err_sst), mur_sst_err, xlab='drifter', ylab='MUR') # scatter plot
data <- data.frame( x=10 * log10(err_sst), y=mur_sst_err) # 2-d histogram
ggplot(data, aes(x=10 * log10(err_sst), y=mur_sst_err) ) +
  geom_bin2d(bins = 200) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()

# Q4: Assuming independence, what is the standard error of the differences? Consult your lecturer if unsure!
# are these standard errors broadly consistent with the observed standard deviations of the differences?
# what might you do next to establish if the uncertainty estimates and standard errors agree?
err_com <- mean(sqrt(err_sst^2 + mur_sst_err^2)) # formula for standard deviation of differences
# next steps, you could do a chi-square test for the variance for example, see here:
# https://www.itl.nist.gov/div898/handbook/eda/section3/eda358.htm#:~:text=A%20chi%2Dsquare%20test%20(%20Snedecor,or%20a%20one%2Dsided%20test.