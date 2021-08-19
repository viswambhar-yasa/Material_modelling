#AUTHOR: VISWAMBHAR REDDY YASA
#MATRICULATION NUMBER :65074
#STOCHASTIC FOR MATERIAL SCIENCE : PROGRAMMING ASSIGNMENT
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc()
#install.packages("stringr")
require(stringr)
require(bitops)
source('./libs.r')
### Installation of important libraries
source('./fcts.r')
### contain all important function required for the programming task
pixel_len=0.001
cat("TASK 03: 3D structure consisting of spherical grains represented using 2D sections \n")
file_path1=file.path('./image3_S04.png')
image3=readImage(file_path1)
display(image3,method = 'raster')

cat("\nTASK 03 (a): Emperical characteristics of the random set \n")
charactertics=estALX(image3,spacing = pixel_len)
#applying transformation to represent 3D CHARACTERISTICS
VF=charactertics[1]
SF=(4/pi)*charactertics[2]
cat('Volume fraction : ')
cat(VF)
cat('\nSpecific surface area : ')
cat(SF)


cat("\n \nTASK 03 (b): segmentation algorithm \n")
?bwlabel
img3_label=bwlabel(image3)
## extracting features from the image
cat('Extracting features \n')
?computeFeatures.moment
img_features = computeFeatures.moment(img3_label)
dim(img_features)
attributes(img_features)$dimnames[[2]]

# considering features which are within observation window
img_features_obser=img_features[img_features[,1]>30 & img_features[,1]<1400 & img_features[,2]>30 & img_features[,2]<1400, 1:3]

disc_centers = cbind(img_features_obser[,1],img_features_obser[,2])/1000
disc_diameters = img_features_obser[,3]/1000 # converting to millimeters 
number_of_grains=length(disc_diameters)

cat('Number of discs within observation window :')
cat(number_of_grains)

cat('\nRange of the disc diameter :[')
range_grain_diameter= range(disc_diameters)
cat(range_grain_diameter)
mean_dia = mean(disc_diameters)

cat(']\nMean of the disc diameter :')
cat(mean_dia)

cat('\nVariance of the disc diameter : ')
sd_dia = sd(disc_diameters)
cat(sd_dia)

cat('\nPlotting histogram of disc diameters')
png(file="03_histogram_disc_dia.png")
disc_count_fr=hist(disc_diameters, breaks=seq(0, 0.2, by=0.01),col=rgb(0.1,0.1,0.7,0.5))$density
dev.off()


cat("\n \nTASK 03 (c): Scheil-Schwartz-Saltykov method\n")

Area_reduced = prod(dim(image3)-2*30)/1000^2 
disc_count_3d=saltykov(Delta=0.001, disc_count_fr, nEM=25)
disc_count_3d
bins=seq(0.01, 0.20, by=0.01)
bins
cat("Plotting histogram for 2d and 3d diameter distributions")
png(file="03_histogram_disc_dia_2d_3d.png")
barplot(rbind(disc_count_fr/sum(disc_count_fr),disc_count_3d/sum(disc_count_3d)),xlab='radius [mm]',ylab='frequency',main="Distribution of diameters in 2d and 3d" ,beside=TRUE, legend.text = c("2D", "3D"),col = c("red", "Blue"),
        border = c("red","Blue") ,args.legend = list(x = "topright"),name=bins)
dev.off()

cat("\n \nTASK 03 (d): estimates of the intensity and mean diameter\n")

intensity = sum(1/disc_diameters)*2/pi/Area_reduced 
cat("Intensity(mean number of balls per unit volume):")
cat(intensity)
cat("\nMean ball diameter 3D:")
mean_diameter = disc_count_fr*pi/2/sum(1/disc_diameters) 
cat(sum(mean_diameter))
