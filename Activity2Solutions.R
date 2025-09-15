# ACTIVITY 2: ANALYSING DRIFTER TRAJECTORIES AND THEIR POWER SPECTRAL DENSITIES

# The goal of this activity is to analyse two drifter trajectories. Learning Outcomes:
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
plot(lon,lat,type = "l",xlab = "longitude", ylab = "latitude") # spaghetti plot

# Q2: Consider 3 segments of length 120 (10 days) starting at points 3301, 3901 and 4501 respectively
# Plot these portions in a different colour on top of your last plot
# Zoom in on longitude range (-72,-65) and latitude rang (33,38)
W<-120 # segment length
plot(lon,lat,type = "l",xlab = "longitude", ylab = "latitude",xlim=c(-72,-65),ylim =c(33,38))
sp<-3301 # starting point 1
lines(lon[sp:(sp+W-1)],lat[sp:(sp+W-1)],type = "l",col = "red")
sp<-3901 # starting point 2
lines(lon[sp:(sp+W-1)],lat[sp:(sp+W-1)],type = "l",col = "blue")
sp<-4501 # starting point 3
lines(lon[sp:(sp+W-1)],lat[sp:(sp+W-1)],type = "l",col = "green")

# Q3: For each segment compute the power spectral density of the (complex-valued) velocities
# Interpret the plots in each case, what do the peaks tell you?
# Show the location of the Coriolis frequency and frequency zero
# For this code, suggested code is given below which you are welcome to use!
W<-120 # segment length
sp<-3301 # starting point 1
freq <- 12*((-W/2):(W/2-1))/W # frequency in units of cycles per day
taper<-c(dpss(W,1,1)$v) # dpss taper of bandwidth 1 (this improves the approximation of the power spectral density)
spec<-fftshift(abs(fft(taper*vel[sp:(sp+W-1)])^2))/12 # power spectral density
plot(freq,10*log10(spec),type = "l",ylab = expression(PSD10log[10] ~ (cm/s)^2~cpd^-1), xlab = "frequency (cycles per day)")
lines(c(0,0),c(min(10*log10(spec))-10,max(10*log10(spec))+10),col="blue")
inert_freq<- coriolis[sp+W/2] # Coriolis frequency at middle of segment
lines(c(inert_freq,inert_freq),c(min(10*log10(spec))-10,max(10*log10(spec))+10),col="red")
legend(2, 30, legend=c("PSD", "Zero frequency","Coriolis frequency"),col=c("black","blue", "red"), lty=1, cex=0.8)

# Q4: Now plot a spectrogram of the complex-valued velocities
# Initially use a window of width 120 data points as in Q2 and Q3
# Then experiment with changing this to see if you can observe the Heisenberg uncertainty principle!
W<-120 # segment length
taper<-c(dpss(W,1,1)$v) # dpss taper of bandwidth 1 (this improves the approximation of the power spectral density)
N<-length(vel) # time series length
STFT<-matrix(0,nrow=N+1-W,ncol=W) # create a matrix for the spectrogram
for (ii in 1:(N+1-W)) {
  STFT[ii, ] <- fftshift(abs(fft(taper*vel[ii:(ii+W-1)]))^2)} # loop over time to populate the matrix with power spectral densities
days <- (W/2+(1:(N+1-W)))/12 # number of days
image_data <- 10*log10(STFT) # plot the STFT on a decibel scale
image(days,freq,image_data,col=rainbow(50))
lines(c(min(days),max(days)),c(0,0),col="black") # add lines for 0 and Coriolis
lines(days,coriolis[(1+W/2):(W/2+N+1-W)],col="green")
legend(300, 5, legend=c("Zero frequency","Coriolis frequency"),col=c("black", "green"), lty=1, cex=0.8)

# Q5: Now try and plot the whole trajectory on a map of Earth - one possibility is to use ggplot2 features
# Try zooming in the plot spanning longitudes (-120,-30) and latitudes (10,60)
# Load some libraries
library(ggplot2)
library(sf)
# Create a data frame of your trajectory
drifter <- data.frame(latitude = lat,longitude = lon)
# Convert your drifter data into a simple features (sf) object
drifter_sf <- st_as_sf(drifter, coords = c("longitude", "latitude"), crs = 4326)
# Plotting the trajectory
ggplot() +
  borders("world", colour = "gray85", fill = "gray80") + # Add the world map
  geom_sf(data = drifter_sf, aes(geometry = geometry), color = "red", size = .01) + # Plot drifter trajectory
  coord_sf(crs = st_crs(4326), xlim = c(-120,-30), ylim = c(10,60)) + # Set the coordinate reference system
  labs(title = "Drifter Trajectory", x = "Longitude", y = "Latitude") +
  theme_minimal()

# Q6: Repeat the analysis on drifterulysses! Can you spot tides in this dataset?!
# This is a length 19704 trajectory of positions sampled every 2 hours starting from 02/01/2005 at 00:07:12 UTC

# Q7: Download the entire dataset (drifters.mat), and create visualisations of all the trajectories
# Try make a plot of average speeds or of average spectra as in the slides, or of some other statistic of your choice?