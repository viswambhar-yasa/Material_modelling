#AUTHOR: VISWAMBHAR REDDY YASA
#MATRICULATION NUMBER :65074
#STOCHASTIC FOR MATERIAL SCIENCE : PROGRAMMING ASSIGNMENT
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc()
source('libs.r')
source('fcts.r')
#install.packages("stringr")
require(stringr)
cat("TASK04: Monte Carlo approach to model Matérn III hard-disc \n")
lamdba=c(4,8,12)
R = 0.06
spacing = 0.01
w=10
n=1
cat('Number of random experimeters :')
cat(n)
for (j in lamdba){
  cat('\n  Simulation is running')
  characteristics_HB3=c(0,0,0)
  for (i in 1:n){
    HD3 = rM3.disc.const(lambda=j, R=R, W=owin(c(0,w), c(0,w)))
    BW.HB3 = digitizeDiscSys(HD3, spacing=spacing) # resolution 'spacing=0.01'
    characteristics_HB3= cbind(characteristics_HB3,estALX(BW.HB3, spacing=spacing))
  }
  cat("\nlamdba :")
  cat(j)
  cat("\n Emperical characteristics of matern hard ball model 3 ")
  
  cat("\n Mean of the area fraction : ")
  Area_mean = mean(characteristics_HB3[1,2:n])
  cat(Area_mean)
  cat("\n variance of the area fraction : ")
  Area_sd = sd(characteristics_HB3[1,2:n])
  cat(Area_sd)
  
  cat("\n Mean of the perimeter fraction : ")
  perimeter_mean = mean(characteristics_HB3[2,2:n])
  cat(perimeter_mean)
  cat("\n variance of the perimeter fraction : ")
  perimeter_sd = sd(characteristics_HB3[2,2:n])
  cat(perimeter_sd)
  
  cat("\n Mean of the Euler Number: ")
  Euler_no_mean = mean(characteristics_HB3[3,2:n])
  cat(Euler_no_mean)
  cat("\n Variance of the Euler Number: ")
  Euler_no_sd = sd(characteristics_HB3[3,2:n])
  cat(Euler_no_sd)
  
  cat("\n Plotting boxplots of the emperical characteristics\n")
  
  file_name1=str_replace_all((paste('./04_area_fraction',trimws(as.character(j)),'.png'))," ","")
  png(file=file_name1)
  title1=paste("Distribution of area fraction for lamda",trimws(as.character(j)))
  boxplot(characteristics_HB3[1,2:n],main=title1,ylab='Area Fraction A[A]')
  dev.off()
  
  file_name2=str_replace_all((paste('./04_perimeter_fraction',trimws(as.character(j)),'.png'))," ","")
  png(file=file_name2)
  title2=paste("Distribution of perimeter fraction for lamda",trimws(as.character(j)))
  boxplot(characteristics_HB3[2,2:n],main=title2,ylab='perimeter L[A]')
  dev.off()
  
  file_name3=str_replace_all((paste('./04_Euler_number',trimws(as.character(j)),'.png'))," ","")
  png(file=file_name3)
  title3=paste("Distribution of Euler Number for lamda",trimws(as.character(j)))
  boxplot(characteristics_HB3[3,2:n],main=title3,ylab='Euler Number X[A]')
  dev.off()
}
