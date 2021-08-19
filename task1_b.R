# AUTHOR : YASA VISWAMBHAR REDDY
# MATRICULATION No: 65074
# STOCHASTICS FOR MATERIAL SCIENCE : PROGRAMMING ASSIGNMENT
cat('PROGRAM IS EXECUTING')
#rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc()
source('./libs.r')
### Installation of important libraries 

cat("TASK 01: Evaluation of Minkowski functions \n")
file_path1=file.path('./image1_S04.png')
file_path2=file.path('./image2_S04.png')
cat("Importing data into R environment \n")
image01<-readImage(file_path1,type="png")
image01@.Data
display(image01,method = 'raster')
image02<-readImage(file_path2,type="png")
image02@.Data
display(image02,method = 'raster')

source('./fcts.r')
### contain all important function required for the programming task
m=0.025 #pixel length
estALX(image01@.Data,spacing=m)
estALX(image02@.Data,spacing=m)
#'estALXFct' is used to evaluate the Minkowski functions
#'estALXFct requires image, m (size of radius), spacing (size of the pixel), ms (minus sampling), shape='disc' (shape of the random object)

cat('Estimation of emperical characteristics like Area, perimeter and Euler Number \n')

image01_char=estALXFct(image01,m=30,spacing=m,shape='disc')
image01_char


image02_char=estALXFct(image02,m=30,spacing=m,shape='disc')
image02_char

cat('Plotting and saving the Minkowski functions  \n')
# show.type : 0 for all , 1 for Area, 2= Area, 3-Euler Number
png(file = "./image01_all.png")
plotALXFct(BW=image01, ALX=image01_char, ALX2=NULL, show.type=0)
dev.off()

png(file = "./image02_all.png")
plotALXFct(BW=image02, ALX=image02_char, ALX2=NULL, show.type=0)
dev.off()

png(file='./01_Compartion_area.png')
plotALXFct(BW=image01, ALX=image01_char,ALX2=image02_char, show.type=1)
legend(0.6,0.8, legend=c("image01", "image02"),col=c("red", "blue"), lty=1:1, cex=1,box.lty=0)
dev.off()

png(file='./01_Compartion_surface.png')
plotALXFct(BW=image01, ALX=image01_char,ALX2=image02_char, show.type=2)
legend(0.6,0.8, legend=c("image01", "image02"),col=c("red", "blue"), lty=1:1, cex=1,box.lty=0)
dev.off()


png(file='./01_Compartion_euler_number.png')
plotALXFct(BW=image01, ALX=image01_char,ALX2=image02_char, show.type=2)
legend(0.6,0.8, legend=c("image01", "image02"),col=c("red", "blue"), lty=1:1, cex=1,box.lty=0)
dev.off()

cat('Plotting of the characteristics is successfull \n')
cat('Execution completed')