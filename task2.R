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

cat("TASK 02: morphological operations and Boolean model \n")
file_path1=file.path('./image1_S04.png')
image1=readImage(file_path1)

cat("Performing morphological opening operation for a list of radius \n")
r_size=c(1,3,5,7,9,11)
cat(r_size)
for(i in r_size){
  B = makeBrush(size=i,shape='disc')
  file_name=str_replace_all((paste('./02_morphological_opening',trimws(as.character(i)),'.png'))," ","")
  png(file=file_name)
  
  morphological_oper=opening(image1, B)
  display(morphological_oper, method="raster", interpolate=FALSE)
  dev.off()
}

n=999
m1=40
sp=0.025
alpha=0.05
cat("\nBuilding Boolean Model for person A \n")
lamda1=4.9
a1=2.5 # length of the rectangle
b1=0.15 # breath of the rectangle
w1=10 # observation window dimensions

#'rBM.rect.const' is obtained from fcts.r  which generates centers and direction angles of rectangles according to a Boolean model
# with intensity 'lambda' and typical grain equal to a uniformly oriented constant rectangle of side lengths 'a' and 'b'
# inside the window 'W' employing plus-sampling correction
XYABP=rBM.rect.const(lamda1,a1,b1,owin(c(0,w1),c(0,w1)))
# transforms a rectangle system 'A XYABP (centre coordinates, side lengths, angle) into a 2D binary image
# with pixel length 'spacing', the rectangular cutout the image represents is given by 'xrange' and 'yrange' 
BW.XYABP_1=digitizeRectSys(XYABP,sp)
model_1=estALXFct(BW=BW.XYABP_1, m=m1, spacing=sp, ms = TRUE)
display(BW.XYABP_1, method="raster")

# Below function which performs simulation of boolean model for the given parameters

simBM.rect.const = function(nrep = 39, lambda, a, b, W = owin(c(0,1),c(0,1)), m=0, spacing=1) {
  sapply(1:nrep, function(k) {
    XYABP = rBM.rect.const(lambda=lambda, a=a, b=b, W=W)
    BW.XYABP = digitizeRectSys(XYABP, spacing=spacing) 
    estALXFct(BW=BW.XYABP, m=m, spacing=spacing, ms = TRUE)[,2]
  })
}
boolean_model_1 = simBM.rect.const(nrep=n, lambda=lamda1, a=a1, b=b1 , W = owin(c(0,w1),c(0,w1)), m=m1, spacing=sp) 
boolean_model_1

cat('Checking the goodness of the fit using envolpe test \n')
fit_test1=globalTest(model_1,boolean_model_1)
cat('Performance for parameter of person A :')
cat(fit_test1)
envolpe_testing1<-globalEnvelopes(boolean_model_1, alpha)

cat('\nPlotting characteristics graph for Person A parameters\n')
png(file='./02_Boolean_model_personA.png')
plot(model_1[,1], model_1[,2], type="l", xlab="r", ylab=expression(paste(A[A],"(r)",sep="")), lwd=2, ylim=c(0.4,1),col='green')
lines(model_1[,1],boolean_model_1[,1],lwd=2, col=2)
lines(model_1[,1],boolean_model_1[,n],lwd=2, col=2)
legend(0.2,0.85, legend=c("upper limit", "lower limit",'sample fit'),col=c("red", "red","green"), lty=1:1, cex=1,box.lty=0)
dev.off()


cat("\nBuilding Boolean Model for person B \n")
lamda2=4.8
a2=1.8 # length of the rectangle
b2=0.13 # breath of the rectangle
w2=10 # observation window dimensions

XYABP_2=rBM.rect.const(lamda2,a2,b2,owin(c(0,w2),c(0,w2)))
# transforms a rectangle system 'A XYABP (centre coordinates, side lengths, angle) into a 2D binary image
# with pixel length 'spacing', the rectangular cutout the image represents is given by 'xrange' and 'yrange' 
BW.XYABP_2=digitizeRectSys(XYABP_2,sp)
model_2=estALXFct(BW=BW.XYABP_2, m=m, spacing=sp, ms = TRUE)
display(BW.XYABP_2, method="raster")

boolean_model_2= simBM.rect.const(nrep=n, lambda=lamda2, a=a2, b=b2 , W = owin(c(0,w2),c(0,w2)), m=m, spacing=sp) 
boolean_model_2

cat('Checking the goodness of the fit using envolpe test \n')
fit_test2=globalTest(model_2,boolean_model_2)
cat('Performance for parameter of person B :')
cat(fit_test2)
envolpe_testing2=globalEnvelopes(boolean_model_2, alpha)

cat('\n Plotting characteristics graph for Person B parameters \n')
png(file='./02_Boolean_model_personB.png')
plot(model_2[,1], model_2[,2], type="l", xlab="r", ylab=expression(paste(A[A],"(r)",sep="")), lwd=2, ylim=c(0.4,1))
lines(model_2[,1],boolean_model_2[,1],lwd=2, col=2)
lines(model_2[,1],boolean_model_2[,n],lwd=2, col=2)
legend(0.2,0.85, legend=c("upper limit", "lower limit",'sample fit'),col=c("red", "red","black"), lty=1:1, cex=1,box.lty=0)
dev.off()