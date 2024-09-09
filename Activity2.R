# ACTIVITY 2: ANALYSING DRIFTER TRAJECTORIES AND THEIR POWER SPECTRAL DENSITIES

# The goal of this activity is to analyse individual drifter trajectories. Learning Outcomes:
# After this computer lab session you should be able to
# - plot trajectories in R, possibly on Earth-based maps
# - calculate and interpret power spectral densities
# - calculate and interpret spectrograms

# Load some libraries (uncomment the installation lines if you have not installed these packages before)
# install.packages("R.matlab")
library(R.matlab)
# install.packages("multitaper")
library(multitaper)
# install.packages("gsignal")
library(gsignal)

# Load data, which first corresponds to a length 6154 trajectory of positions sampled every 2 hours starting from 13/10/2024 at 02:35:31 UTC
data <- readMat('drifterbetty.mat')
betty = data$drifterbetty[, , 1]
lat = betty$lat # drifter latitudes
lon = betty$lon # drifter longitudes
vel = betty$cv # drifter complex-valued velocities
coriolis = -sign(lat)*24*betty$f # drifter Coriolis frequencies in cycles per day

# Q1: Using the plot function in R, Plot the entire trajectory as a line or strand of "spaghetti" - what features do you see?


# Q2: Consider 2 segments of length 120 (10 days) starting at points 3301, 3901 and 4501 respectively
# Plot these portions in a different colour on top of your last plot
# Zoom in on longitude range (-72,-65) and latitude rang (33,38)


# Q3: For each segment compute the power spectral density of the (complex-valued) velocities
# Interpret the plots in each case, what do the peaks tell you?
# Show the location of the Coriolis frequency and frequency zero
# For this code, suggested code is given below which you are welcome to use!
W<-120
sp<-3301
freq <- 12*((-W/2):(W/2-1))/W # frequency in units of cycles per day
taper<-c(dpss(W,1,1)$v) # dpss taper of bandwidth 1
spec<-fftshift(abs(fft(taper*vel[sp:(sp+W-1)])^2))/12 # power spectral density
plot(freq,10*log10(spec),type = "l",ylab = expression(PSD10log[10] ~ (cm/s)^2~cpd^-1), xlab = "frequency (cycles per day)")
lines(c(0,0),c(min(10*log10(spec))-10,max(10*log10(spec))+10),col="blue")
inert_freq<- coriolis[sp+W/2] # Coriolis frequency at middle of segment
lines(c(inert_freq,inert_freq),c(min(10*log10(spec))-10,max(10*log10(spec))+10),col="red")
legend(2, 30, legend=c("PSD", "Zero frequency","Coriolis frequency"),col=c("black","blue", "red"), lty=1, cex=0.8)

# Q4: Now plot a spectrogram of the complex-valued velocities
# Initially use a window of width 120 data points as in Q2 and Q3
# Then experiment with changing this to see if you can observe the Heisenberg uncertainty principle!


# Q5: Now try and plot the whole trajectory on a map of Earth - one possibility is to use ggplot2 features
# Try zooming in the plot spanning longitudes (-120,-30) and latitudes (10,60)


# Q6: Repeat the analysis on drifterulysses! Can you spot tides in this dataset?!
# This is a length 19704 trajectory of positions sampled every 2 hours starting from 02/01/2005 at 00:07:12 UTC
# You can also try drifterinti, drifterisis, driftermagellan, drifternansen, all are quite different!


# Q7: (Optional!) Download the entire dataset (ask your lecturer for the link!), and create visualisations of all the trajectories
# Try make a plot of average speeds or of average spectra as in the slides, or of some other statistic of your choice?