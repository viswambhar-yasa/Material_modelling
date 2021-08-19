########
# fcts.r
# useful R functions
#
# You can make R to know all subsequent functions by running source('fcts.r').
########

###
# Estimation of quermass densities and Minkowski functions
###

estALX<-function(B, spacing=1) {
  # estimates A_A, L_A and chi_A from a 2D binary image 'B' with pixel length 'spacing' 
  require(bitops, quietly=TRUE)
  acot<-function(x) pi/2-atan(x)
  F1<-matrix(c(1,4,2,8),2,2,byrow=T)
  s<-dim(B)
  n0<-prod(s-1)
  G1<-B[1:(s[1]-1),1:(s[2]-1)]*F1[1,1]+B[2:s[1],1:(s[2]-1)]*F1[2,1]+B[1:(s[1]-1),2:s[2]]*F1[1,2]+B[2:s[1],2:s[2]]*F1[2,2]
  
  h<-rep(0,16)
  for(k in 0:15) h[k+1]<-sum(G1==k)
  
  theta<-c(0, acot(1), pi/2, pi-acot(1), pi, pi+acot(1), 3*pi/2, 2*pi-acot(1))
  r<-c(1, sqrt(2), 1,  sqrt(2), 1, sqrt(2), 1, sqrt(2))
  kappa0<-c(1, 1, 4, 4, 8, 8, 2, 2)
  kappa1<-c(4, 8, 8, 2, 2, 1, 1, 4)
  k<-1:16
  pA<-rep(0,8)
  for(l in 1:8){
    pA[l]<-sum((bitOr(k-1,kappa0[l])==(k-1))*(bitAnd(k-1,kappa1[l])==0)*h)/r[l]/n0
  }
  
  cc<-rep(0,8)
  cc[1]=(2*pi-theta[8]+theta[2])/4/pi
  cc[8]=(2*pi-theta[7]+theta[1])/4/pi
  for(l in 2:7) cc[l]=(theta[l+1]-theta[l-1])/4/pi
  
  AA<-sum((bitOr(k-1,1)==(k-1))*h)/n0;
  LA<-pi*sum(pA*cc)
  chiA<-(h[2]+h[3]+h[5]+h[9]-h[8]-h[12]-h[14]-h[15])/4/n0
  
  return(c(AA=AA,LA=LA,XA=chiA)/spacing^c(0,1,2))
}

# The following function 'estALXFct' performs the estimation of the Minkowski functions A_A(r), L_A(r) and chi_A(r)
# for the black-white image 'BW' for arguments r from (0,1,...,m)*spacing.
# 'spacing' is the real pixel length.
# By 'ms' minus-sampling correction can be employed ('ms=TRUE') or supressed ('ms=FALSE').
estALXFct <- function(BW, m=0, spacing=1, ms=FALSE, shape='disc') {
  require(EBImage, quietly=TRUE)
  rr <- 0:m
  d <- dim(BW)
  cbind(r=rr*spacing, t(sapply(rr, function(r) estALX(dilate(BW, makeBrush(size=2*r+1, shape=shape))[(1+ms*r):(d[1]-ms*r),(1+ms*r):(d[2]-ms*r)]@.Data, spacing=spacing))))
}

plotALXFct <- function(BW, ALX, ALX2=NULL, show.type=0) {
  oldpar <- par(no.readonly = TRUE)
  if(show.type == 0) {
    layout(matrix(1:4,2,byrow=TRUE))
    display(BW, method="raster", margin=rep(10, 4))
  }
  par(mar=c(4.2,4.2,2,2))
  #
  if(show.type %in% c(0,1)) {
    plot(ALX[,1], ALX[,2], type="n", main=expression(paste("1st Minkowski function ",A[A],"(r)",sep="")), xlab="r", ylab=expression(paste(A[A],"(r)",sep="")))
    abline(h=1, col=8)
    lines(ALX[,1], ALX[,2], lwd=2, col=2)
    if(!is.null(ALX2)) lines(ALX2[,1], ALX2[,2], lwd=2, col=4)
  }
  #
  if(show.type %in% c(0,2)) {
    plot(ALX[,1],ALX[,3],type="n",main=expression(paste("2nd Minkowski function ",L[A],"(r)",sep="")),xlab="r",ylab=expression(paste(L[A],"(r)",sep="")),lwd=2,col=2)
    abline(h=0, col=8)
    lines(ALX[,1], ALX[,3], lwd=2, col=2)
    if(!is.null(ALX2)) lines(ALX2[,1], ALX2[,3], lwd=2, col=4)
  }
  #
  if(show.type %in% c(0,3)) {
    plot(ALX[,1],ALX[,4],type="n",main=expression(paste("3rd Minkowski function ",chi[A],"(r)",sep="")),xlab="r",ylab=expression(paste(chi[A],"(r)",sep="")),lwd=2,col=2)
    abline(h=0, col=8)
    lines(ALX[,1], ALX[,4], lwd=2, col=2)
    if(!is.null(ALX2)) lines(ALX2[,1], ALX2[,4], lwd=2, col=4)
  }
  #
  if(show.type == 0) layout(1)
  par(oldpar)
}

###
# Global rank envlope test and global rank envelopes
###

globalTest<-function(T0, Tmod) { 
  # Determines the global rank envelope test out of dim(Tmod)[2] repetitions 
  # of some functional characteristic.
  # 'Tmod' is a matrix where each column corresponds to one repetition.
  # 'T0' is a vector representing the functional characteristic of the data.
  stopifnot(is.matrix(Tmod))
  nrep <- dim(Tmod)[2]
  TT <- cbind(T0, Tmod)
  RRmin <- t(apply(TT, 1, rank))
  RRmax <- (nrep+2) - RRmin
  RR <- pmin(RRmin, RRmax)
  R <- apply(RR, 2, min)
  plow <- 1 - sum(R[-1]>=R[1]) / (nrep+1)
  pupp <- 1 - sum(R[-1]>R[1]) / (nrep+1)
  return(mean(c(plow,pupp)))
}

globalEnvelopes <- function(Tmod, alpha) 
{
  # Determines the global rank envelope out of dim(Tmod)[2] repetitions 
  # of some functional summary statistic.
  # 'Tmod' is a matrix where each column corresponds to one repetition.
  # 'alpha' is the significance level.
  stopifnot(is.numeric(alpha) && length(alpha) == 1 && alpha > 0 && alpha < 1)
  stopifnot(is.matrix(Tmod))
  nrep <- dim(Tmod)[2]
  RRmin <- t(apply(Tmod, 1, rank))
  RRmax <- (nrep+1) - RRmin
  RR <- pmin(RRmin, RRmax)
  R <- apply(RR, 2, min)
  k <- 1
  while(sum((R<(k+1))) <= alpha*nrep) 
    k <- k+1
  Tlow <- apply(Tmod, 1, sort)[k,]
  Tupp <- apply(Tmod, 1, sort, decreasing=TRUE)[k,]
  return(cbind(Tlow, Tupp))
}

###
# Boolean model with rectangles
###

# The following function 'rBM.rect.const' generates centers and direction angles of rectabgles according to a Boolean model
# with intensity 'lambda' and typical grain equal to a uniformly oriented constant rectangle of side lengths 'a' and 'b'
# inside the window 'W' employing plus-sampling correction.
rBM.rect.const <- function(lambda, a, b, W = owin(c(0,1),c(0,1))) {
  require(spatstat, quietly=TRUE)
  stopifnot("lambda must be positive"= lambda>0)
  stopifnot("a must be positive"= a>0)
  stopifnot("b must be positive"= a>0)
  R <- max(c(a,b))
  W.plus <- owin(W$xrange+c(-R,R), W$yrange+c(-R,R))   # make enlarged window W oplus [-R,R]x[-R,R]
  xy <- rpoispp(lambda=lambda, win=W.plus)             # generate homogeneous Poisson point process
  A <- rep(a, xy$n)                                    # first side length
  B <- rep(b, xy$n)                                    # second side length
  Phi <- runif(xy$n, 0, 2*pi)                          # direction angle
  XYABP <- cbind(X=xy$x, Y=xy$y, A=A, B=B, Phi=Phi)
  attributes(XYABP)$window <- W                        # ensures that also the original window is returned
  return(XYABP)
}

digitizeRectSys <- function(M, spacing, xrange = attributes(XYABP)$window$xrange, yrange = attributes(XYABP)$window$yrange) {
  # transforms a rectangle system 'M' (centre coordinates, side lengths, angle) into a 2D binary image
  # with pixel length 'spacing', the rectangular cutout the image represents is given by 'xrange' and 'yrange'  
  library(spatstat)
  Wdigi <- as.mask(owin(xrange=xrange, yrange=yrange), eps=spacing)
  n <- length(M[,1])
  xCoord <- raster.x(Wdigi)
  yCoord <- raster.y(Wdigi)
  Wdim <- Wdigi[[4]]
  BB <- !Wdigi[[10]]
  for(i in 1:n) {
    imin <- max(c(1, (M[i,1]-0.5*M[i,3]+0.5*spacing)%/%spacing))
    jmin <- max(c(1, (M[i,2]-0.5*M[i,3]+0.5*spacing)%/%spacing))
    imax <- min(c(Wdim[1], (M[i,1]+0.5*M[i,3]+0.5*spacing)%/%spacing))
    jmax <- min(c(Wdim[2], (M[i,2]+0.5*M[i,3]+0.5*spacing)%/%spacing))
    if(imax>0 && jmax>0 && imin<=Wdim[1] && jmin<=Wdim[2]) BB[imin:imax,jmin:jmax] <- pmax(BB[imin:imax,jmin:jmax],t((abs((xCoord[jmin:jmax,imin:imax]-M[i,1])*cos(M[i,5])+(yCoord[jmin:jmax,imin:imax]-M[i,2])*sin(M[i,5]))<=M[i,3]*0.5)*(abs((xCoord[jmin:jmax,imin:imax]-M[i,1])*sin(M[i,5])-(yCoord[jmin:jmax,imin:imax]-M[i,2])*cos(M[i,5]))<=M[i,4]*0.5)))
  }
  return(BB*1.0)
} 

###
# Saltykov's method
###

kernelsp0<-function(u,s) {
  # auxiliary kernel function
  if(s<0.0) return(u) 
  if(s<u) return(sqrt(u*u-s*s)) 
  return(0.0)
}

emalgor<-function(nEM,p,nn,NN) {
  # Given incomplete data 'nn', the matrix 'p', and an initial solution 'NN' 
  # for the parameter vector, this routine performs the EM algorithm
  # where 'nEM' is the number of EM steps.
  k <- length(nn)
  q <- apply(p,2,sum)
  for(count in 1:nEM)
  {
    r <- p%*%NN
    rid <- (1:k)[r>0]
    for(i in 1:k)
    {  
      s <- sum(p[rid,i]*nn[rid]/r[rid])
      if(q[i]>0) NN[i] <- NN[i]*s/q[i];
    }
  }
  return(NN)
}   

saltykov<-function(nn,Delta,nEM=32) {
  # Given the vector 'nn' of absolute frequencies of disc diameters (per unit area)
  # observed in a planar section (planar sampling in a planar section), the
  # vector 'NN' of absolute frequencies of ball diameters (per unit volume) is returned.
  # 'nEM' is the number of iterations in the EM algorithm.
  # 'Delta' is the bin width w.r.t. the diameters.
  k <- length(nn)
  NN <- rep(0,k)
  p <- matrix(0,k,k)
  for(i in 1:k) for(j in 1:k) p[j,i] <- kernelsp0(i,j-1)-kernelsp0(i,j)
  NN <- nn+1e-6
  NN <- emalgor(nEM,p,nn,NN)
  NN <- NN/Delta
  return(NN)
}

###
# Matern III model
###

# Matern III - Perfect Simulation

isInBox <- function(x, box) min((x>=0) & (x<=box))

isInBallRadii <- function(x, z, d) sum((x[1:d]-z[1:d])^2)<=(x[d+1]+z[d+1])^2

isInUnionRadii <- function(x, V, d) max(apply(matrix(V[,1:(d+1)],ncol=d+1),1,isInBallRadii,z=x[1:(d+1)],d=d) & (x[d+2]<=V[,d+2]))

rMatern3Radii<-function(lambda, R, box, p)
{
  X<-c()
  d<-length(box)
  Rmax<-max(R)
  n<-rpois(1,lambda*prod(box))
  Y<-matrix(runif(d*n)*box,ncol=d,byrow=TRUE)
  T<-runif(n)
  V<-c()
  J<-apply((rmultinom(n,1,p)>0),2,order,decreasing=TRUE)[1,] 
  if(n>0) onceagain<-max(apply(Y,1,isInBox,box=box)) else onceagain<-FALSE
  while(onceagain)
  {
    imin<-order(T)[1]
    z<-Y[imin,]
    t<-T[imin]
    j<-J[imin]
    r<-R[j]
    rR<-r+Rmax
    ok<-0
    if(!isInBox(z-rR,box-2*rR)) 
    {
      nn<-rpois(1,(2*rR)^d*lambda*t)
      if(nn>0)
      {
        YY<-matrix(runif(d*nn,-rR,rR)+z,ncol=d,byrow=TRUE)
        TT<-t*runif(nn)
        JJ<-apply((rmultinom(nn,1,p)>0),2,order,decreasing=TRUE)[1,]
        RR<-R[JJ]
        # points outside box and inside B(z,r+R)
        i<-(1:(nn+1))[c((!apply(YY,1,isInBox,box)) & apply(matrix(c(YY,RR),ncol=d+1),1,isInBallRadii,z=c(z,r),d=d),FALSE)]
        # points outside V
        ii<-c()
        if(length(V)>0 & length(i)>0) ii<-c(i,0)[c(!apply(matrix(c(YY[i,],RR[i],TT[i]),ncol=d+2),1,isInUnionRadii,V=V,d=d),FALSE)]
        V<-rbind(V,c(z,r,t))
        if(length(ii)>0)
        {
          Y<-rbind(Y,YY[ii,])
          T<-c(T,TT[ii])
          J<-c(J,JJ[ii])
        }
        else ok<-1 # no further point
      }
      else ok<-1 # no further point
    }
    else ok<-1 # no further point
    if(ok)	
    {
      idx<-(1:length(Y[,1]))[apply((t(Y)-z)^2,2,sum)>(r+R[J])^2]
      Y<-matrix(Y[idx,],ncol=d)
      T<-T[idx]
      J<-J[idx]
      if(isInBox(z,box)) X<-rbind(X,c(z,r))
    }
    if(length(Y[,1])>0) onceagain<-max(apply(Y,1,isInBox,box=box)) else onceagain<-FALSE
  }
  rownames(X)<-NULL
  return(X)
}

# The following function 'rM3.disc.const' generates centers and radii of discs according to a Matern III model
# with intensity 'lambda' and typical grain equal to a constant disc of radius 'R' inside the window 'W'.
rM3.disc.const <- function(lambda, R, W = owin(c(0,1),c(0,1))) {
  R <- rep(R ,2)
  box <- c(W$xrange[2]-W$xrange[1], W$yrange[2]-W$yrange[1])
  XYR <- rMatern3Radii(lambda=lambda, R=R, box=box, p=c(1,0))
  XYR[,1] <- XYR[,1] + W$xrange[1]
  XYR[,2] <- XYR[,2] + W$yrange[1]
  colnames(XYR) <- c("X", "Y", "R")
  attributes(XYR)$window <- W
  return(XYR)
}

digitizeDiscSys <- function(XYR, spacing, xrange = attributes(XYR)$window$xrange, yrange = attributes(XYR)$window$yrange) {
  # transforms a disc system 'XYR' (centre coordinates, radii) into a 2D binary image
  # with pixel length 'spacing', the rectangular cutout the image represents is given by 'xrange' and 'yrange'  
  library(spatstat)
  Wdigi <- as.mask(owin(xrange=xrange, yrange=yrange), eps=spacing)
  n <- length(XYR[,1])
  xCoord <- raster.x(Wdigi)
  yCoord <- raster.y(Wdigi)
  Wdim <- Wdigi[[4]]
  BB <- !Wdigi[[10]]
  for(i in 1:n) {
    imin <- max(c(1, (XYR[i,1]-XYR[i,3]+0.5*spacing)%/%spacing))
    jmin <- max(c(1, (XYR[i,2]-XYR[i,3]+0.5*spacing)%/%spacing))
    imax <- min(c(Wdim[1], (XYR[i,1]+XYR[i,3]+0.5*spacing)%/%spacing))
    jmax <- min(c(Wdim[2], (XYR[i,2]+XYR[i,3]+0.5*spacing)%/%spacing))
    if(imax>0 && jmax>0 && imin<=Wdim[1] && jmin<=Wdim[2]) BB[imin:imax,jmin:jmax] <- pmax(BB[imin:imax,jmin:jmax], t((xCoord[jmin:jmax,imin:imax]-XYR[i,1])^2+(yCoord[jmin:jmax,imin:imax]-XYR[i,2])^2 <= XYR[i,3]^2))
  }
  return(BB*1.0)
} 