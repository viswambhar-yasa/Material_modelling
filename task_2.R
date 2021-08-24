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

sim_rec_boolean_modl = function(nrep = 39, lambda, a, b, W = owin(c(0,1),c(0,1)), m=0, spacing=1) {
  sapply(1:nrep, function(k) {
    XYABP = rBM.rect.const(lambda=lambda, a=a, b=b, W=W)
    BW.XYABP = digitizeRectSys(XYABP, spacing=spacing) 
    estALXFct(BW=BW.XYABP, m=m, spacing=spacing, ms = TRUE)[,2]
  })
}

n=999
m1=30
sp=0.025
alpha=0.05
cat("\nBuilding Boolean Model for person A \n")
lamda1=4.9
a1=2.5 # length of the rectangle
b1=0.15 # breath of the rectangle
w1=10 # observation window dimensions

XYABP = rBM.rect.const(lambda=lamda1, a=a1,b=b1, W=owin(c(0,w1), c(0,w1)))
BW.XYABP = digitizeRectSys(XYABP, spacing=0.025)
display(BW.XYABP, method="raster")
cat('Simulation is running')
boolean_modelA = sim_rec_boolean_modl(nrep=n, lambda=lamda1, a=a1, b=b1 , W = owin(c(0,w1),c(0,w1)), m=m1, spacing=sp) 
cat('Simulation ran successfully')
enevelop_test_A <- globalEnvelopes(boolean_modelA, alpha=0.05)
A_fit = estALXFct(BW=BW.XYABP, m=m1, spacing=sp, ms = TRUE)[,1:2]
cat('Checking the goodness of the fit using envolpe test P= ')
globalTest(AA.data1[,2], boolean_modelA)

cat('\nPlotting characteristics graph for Person A parameters\n')
png(file='./02_Boolean_model_personA.png')
plot(A_fit[,1], A_fit[,2], type="l", xlab="r", ylab=expression(paste(A[A],"(r)",sep="")), lwd=2, ylim=c(0.4,1),xlim=c(0,0.5))
lines(A_fit[,1], enevelop_test_A[,1], lwd=2, col=2)
lines(A_fit[,1], enevelop_test_A[,2], lwd=2, col=2)
legend(0.3,0.85, legend=c("upper limit", "lower limit",'sample fit'),col=c("red", "red","green"), lty=1:1, cex=1,box.lty=0)
dev.off()


cat("\nBuilding Boolean Model for person B \n")
lamda2=4.8
a2=1.8 # length of the rectangle
b2=0.13 # breath of the rectangle
w2=10 # observation window dimensions
m2= 25

XYABP2 = rBM.rect.const(lambda=lamda2, a=a2,b=b2, W=owin(c(0,w2), c(0,w2))) # The rectangular window is [0,10] x [0,10].
BW.XYABP2 = digitizeRectSys(XYABP2, spacing=sp)
display(BW.XYABP2, method="raster")
cat('Simulation is running')
boolean_modelB = sim_rec_boolean_modl(nrep=n, lambda=lamda2, a=a2, b=b2 , W = owin(c(w2),c(0,w2)), m=m2, spacing=sp) 
enevelop_test_B <- globalEnvelopes(boolean_modelB, alpha=0.05)
B_fit = estALXFct(BW=BW.XYABP2, m=m2, spacing=sp, ms = TRUE)[,1:2]
globalTest(B_fit[,2], boolean_modelB)

cat('\nPlotting characteristics graph for Person B parameters\n')
png(file='./02_Boolean_model_personB.png')
plot(B_fit[,1], B_fit[,2], type="l", xlab="r", ylab=expression(paste(A[A],"(r)",sep="")), lwd=2, ylim=c(0.4,1),xlim=c(0,0.5))
lines(B_fit[,1], enevelop_test_B[,1], lwd=2, col=2)
lines(B_fit[,1], enevelop_test_B[,2], lwd=2, col=2)
dev.off()