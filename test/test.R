debug(ncvreg)
library(glmnet)
x=matrix(rnorm(100*20),100,20)
y=rnorm(100)

std <- .Call("standardize", x,F,PACKAGE="bar")
XX <- std[[1]]
center <- std[[2]]
scale <- std[[3]]

ncvreg(x,y,penalty = "SCAD",intercept = T, stand=T,lambda=0.1)
scadreg(x,y,lambdav=0.1)
scadlseq_c(x,y,lambdav=0.1,gam=3.7,1e-5)

scadreg(XX,yy,lambdav=0.1)[[1]]-b
scadlseq_c(XX,yy,lambdav=0.1,gam=3.7,1e-5)
sx <- XX


n  <- 100
p <- 10
x=matrix(rnorm(n*p),n,p)
y=x %*% c(1,1,1,rep(0,p-3))+rnorm(n)
std <- .Call("standardize", x,F,PACKAGE="bar")
XX <- std[[1]]
a1 <- ridreg1_c(XX, y, 0.1, rep(1, 2000), 1e-10)
a2 <- ridreg_c(XX, y, 0.1, rep(1, 2000), 1e-10)
#microbenchmark(ridreg_c(XX, y, 0.1, rep(1, 2000), 1e-10),ridreg1_c(XX, y, 0.1, rep(1, 2000), 1e-10))
A <-t(XX) %*% XX
B <- t(XX) %*% y
system.time(solve(A+diag(rep(1:2000))) %*% B)
dim(XX)
gamm <- diag(1:2000/2000)
#1
system.time(solve(A+diag(rep(1:2000))) %*% B)
system.time(gamm %*% solve(gamm%*%A%*% gamm +diag(rep(0.1,2000)))%*% gamm %*% B)
system.time({dd <- gamm %*% ll$vectors; dd %*% diag(1/(ll$values+1:2000/2000)) %*% t(dd) %*% B})


system.time(ridregM1_c(A,B,0.1,1:p/p,1e-5)[1])
ridregM2_c(A,B,0.1,1:p/p,n,1e-5)[1]
system.time(ridregM3_c(XX,B,0.1,1:p/p,n,1e-5)[1])
system.time(ridregM4_c(XX,B,0.1,1:p/p,min(n,p),1e-5)[1])

#microbenchmark::microbenchmark(ridregM1_c(A,B,0.1,1:p/p,1e-5),ridregM4_c(XX,B,0.1,1:p/p,min(n,p),1e-5))

we <- runif(p)

bin0 <- ginv(A+diag(rep(1,p))) %*% B
bin0 <- ginv(t(x) %*% x+diag(rep(1,p))) %*% t(x) %*% y
eigdecom <- eigs_sym(A,100)

b1 <- ridregM4_c(XX,B,0.1,we,n,1e-5)
b2 <- ridregM1_c(A,B,0.1,we,1e-5)

b3 <- ridregM3_c(XX,B,0.1,we,min(p,n),1e-5)
b30 <- ridregM3_c(t(eigdecom$vectors %*% diag(sqrt(eigdecom$values))),B,0.1,we,n,1e-5)

system.time(t1<- ridreg0_c(XX[,1:200],y,rep(0,200),0.1,we,1e-10))


b11 <- ridregi_c(XX,y,0.1,bin0,1e-7)

b12 <- barrep_c(XX,y,0.1,bin0,1e-7,1)

aset1 <- which(abs(b11)>1e-10)
solve(A[aset1,aset1]+0.1*diag(1/b11[aset1]^2),B[aset1])-b11[aset1]
sum((solve(A[aset1,aset1]+0.1*diag(1/b11[aset1]^2),B[aset1])-b11[aset1])^2)
sum(abs(solve(A[aset1,aset1]+0.1*diag(1/b11[aset1]^2),B[aset1])-b11[aset1]))
sqrt(sum((updatefun_c(XX,y,0.1,b11)-b11)^2))
sum((b11-c(1,1,1,rep(0,1000-3)))^2)

aset2 <- which(abs(b12)>1e-10)
solve(A[aset2,aset2]+0.1*diag(1/b12[aset2]^2),B[aset2])-b12[aset2]
sum((solve(A[aset2,aset2]+0.1*diag(1/b12[aset2]^2),B[aset2])-b12[aset2])^2)
sum(abs(solve(A[aset2,aset2]+0.1*diag(1/b12[aset2]^2),B[aset2])-b12[aset2]))
sqrt(sum((updatefun_c(XX,y,0.1,b12)-b12)^2))
sum((b12-c(1,1,1,rep(0,1000-3)))^2)

microbenchmark(ridregi_c(XX,y,0.1,bin0,1e-7))
microbenchmark(ridregi_c(XX,y,0.1,bin0,1e-7), barrep_c(XX,y,0.1,bin0,1e-7))
b11

barrep_c(XX,y,0.1,bin0,1e-7)

b0 <- barrep_c(x,y,0.1,bin0,1e-7,0)
b1 <- barrep_c(x,y,0.1,bin0,1e-7,1)
b2 <- itercd_c(x,y,p,n,0.1,bin0,1e-7)

b0 <- barrep_c(XX,y,0.1,bin0,1e-7,0)
b1 <- barrep_c(XX,y,0.1,bin0,1e-7,1)
b2 <- itercd_c(XX,y,p,n,0.1,bin0,1e-7)
microbenchmark(b0 <- barrep_c(XX,y,0.1,bin0,1e-7,0),
b1 <- barrep_c(XX,y,0.1,bin0,1e-7,1),
b2 <- itercd_c(XX,y,p,n,0.1,bin0,1e-7))
sum((b2-c(1,1,1,rep(0,p-3)))^2)
sum((b1-c(1,1,1,rep(0,p-3)))^2)
sqrt(sum((updatefun_c(XX,y,0.1,b1)-b1)^2))
sqrt(sum((updatefun_c(XX,y,0.1,b2)-b2)^2))
which(abs(b1)>1e-20)

foo <- function(n,p,k,l){
  b <- c(rep(1,k),rep(0,p-k))
  x <- matrix(rnorm(n*p),n,p)
  y <- x %*% b+rnorm(n)
  bin0 <- ginv(t(x) %*% x+diag(rep(1,p))) %*% t(x) %*% y
  lmax <- max((t(x) %*% y)^2 / diag(t(x) %*% x))
  lseq <- exp(seq(log(lmax),log(lmax*0.0001),length.out=50))
  res <- matrix(0,nrow=50,ncol=8)
  for( i in 1:50){
      barrep_c(x,y,lseq[i],bin0,1e-7)
    bic0 <- n*log(sum((y-x%*%b0)^2)/n)+sum(b0!=0)*log(n)
    b2 <- itercd_c(x,y,p,n,lseq[i],bin0,1e-8)
    bic2 <- n*log(sum((y-x%*%b2)^2)/n)+sum(b2!=0)*log(n)
    #browser()
    res[i,] <- c(sum(b0[1:k]!=0),sum(b0[-(1:k)]!=0),sum((b0-b)^2),bic0,sum(b2[1:k]!=0),sum(b2[-(1:k)]!=0),sum((b2-b)^2),bic2)
  }
  browser()
  res1 <- c(res[which(res[,4]==min(res[,4]))[1],1:4],res[which(res[,8]==min(res[,8]))[1],5:8])
  return(res1)
}



library(doParallel)
stopCluster(cl)
cl <- makeCluster(6)
registerDoParallel(cl)

rr <- foreach(i=1:100,.packages = c("MASS","bar"),.combine="rbind") %dopar% foo(10,100,3,0.1)

rr <- foreach(i=1:500,.packages = c("MASS","bar"),.combine="rbind") %dopar% foo(100,10,3,0.1)

#rr <- replicate(100,foo(100,10,3,0.1))
colMeans(rr)

foo(100,10,3,0.1)

tr <- foreach(i=exp(seq(log(10),log(10*0.0001),length.out=50)),.packages = c("MASS","bar"),.combine="cbind") %dopar% {barrep_c(x,y,i,bin0,1e-7,0)}
barrep_c(x,y,0.001,bin0,1e-7,0)

foo1 <- function(n,p,k,l){
  b <- c(rep(1,k),rep(0,p-k))
  x <- matrix(rnorm(n*p),n,p)
  y <- x %*% b+t(rmvnorm(1,rep(0,n),diag(rep(1,n))))*3
  bin0 <- ginv(t(x) %*% x+diag(rep(1,p))) %*% t(x) %*% y
  lmax <- max((t(x) %*% y)^2 / diag(t(x) %*% x))
  lseq <- exp(seq(log(lmax),log(lmax*0.01),length.out=50))
  #res <- matrix(0,nrow=50,ncol=8)
  #browser()
  cv0 <- cv_test1_c(x,y,lseq,bin0,1e-8,5)
  cv2 <- cv_test2_c(x,y,lseq,bin0,1e-9,5)
  b0 <- cv0[[1]]
  b2 <- cv2[[1]]
  browser()
  res1 <- c(sum(abs(b0[1:k])>1e-10),sum(abs(b0[-(1:k)])>1e-10),sum((b0-b)^2),sum(abs(b2[1:k])>1e-10),sum(abs(b2[-(1:k)])>1e-10),sum((b2-b)^2))
  return(res1)
}
foo1(50,100,3,0.1)

foo1(30,100,3,0.1)

library(doParallel)
stopCluster(cl)
cl <- makeCluster(7)
registerDoParallel(cl)
rr <- foreach(i=1:500,.packages = c("MASS","bar"),.combine="rbind") %dopar% foo1(50,10,3,0.1)
rr <- foreach(i=1:500,.packages = c("MASS","bar"),.combine="rbind") %dopar% foo1(30,50,3,0.1)
rr <- foreach(i=1:100,.packages = c("MASS","bar","mvtnorm"),.combine="rbind") %dopar% foo1(30,100,3,0.1)
cv_test1_c(x,y,seq(1,0.01,length.out = 10),bin0,1e-7,5)


n <- 30
p <- 100
x <- matrix(rnorm(n*p),n,p)
y <- x %*% c(1,1,1,rep(0,p-3))+rnorm(n)
bin0 <- ginv(t(x) %*% x+diag(rep(1,p))) %*% t(x) %*% y
tt <- cv_test3_c(x,y,exp(seq(log(10),log(0.001),length.out = 50)),bin0,1e-7,5)
t(tt[[2]]-tt[[5]])
t(tt[[3]])
t(tt[[6]])

as.vector(itercd_c(x,y,p,n,0.2442053,bin0,1e-8))
as.vector(barrep_c(x,y,0.2442053,bin0,1e-8))


tt <- cv_test3_c(x,y,exp(seq(log(0.001),log(10),length.out = 10)),exp(seq(log(100),log(1),length.out=10)),bin0,1e-7,5)



foo2 <- function(n,p,k,l){
  b <- c(rep(1,k),rep(0,p-k))
  x <- rmvnorm(n,rep(0,p),toeplitz(0.3^(1:p)))
  y <- x %*% b+t(rmvnorm(1,rep(0,n),diag(1,n)))*1
  A <- t(x) %*% x

  bin0 <- ginv(t(x) %*% x+diag(rep(1,p))) %*% t(x) %*% y
  lmax <- max((t(x) %*% y)^2 / diag(t(x) %*% x))
  lseq <- exp(seq(log(lmax*0.001),log(lmax),length.out=50))
  
  xiv <- exp(seq(log(100),log(1),length.out=10))
  #res <- matrix(0,nrow=50,ncol=8)
  #browser()
  cvres <- cv_test3_c(x,y,lseq,xiv,bin0,1e-7,5)
  #cv0 <- cv_test1_c(x,y,lseq,bin0,1e-8,5)
  #cv2 <- cv_test2_c(x,y,lseq,bin0,1e-9,5)
  b0 <- cvres[[1]]
  b2 <- cvres[[4]]
  browser()
  res1 <- c(sum(abs(b0[1:k])>1e-10),sum(abs(b0[-(1:k)])>1e-10),sum((b0-b)^2),sum(abs(b2[1:k])>1e-10),sum(abs(b2[-(1:k)])>1e-10),sum((b2-b)^2))
  return(res1)
}

foo3 <- function(n,p,k,l){
  b <- c(rep(1,k),rep(0,p-k))
  x <- matrix(rnorm(n*p),n,p)
  y <- x %*% b+t(rmvnorm(1,rep(0,n),diag(rep(1,n))))*1
  A <- t(x) %*% x
  eidecom <- eigs_sym(A,min(n,p))
  xdec <- diag(sqrt(eidecom$values)) %*% t(eidecom$vectors)
  ydec <- t(x) %*% y
  bin0 <- ginv(t(x) %*% x+diag(rep(1,p))) %*% t(x) %*% y
  lmax <- max((t(x) %*% y)^2 / diag(t(x) %*% x))
  lseq <- exp(seq(log(lmax*0.001),log(lmax),length.out=50))
  
  xiv <- exp(seq(log(100),log(1),length.out=10))
  #res <- matrix(0,nrow=50,ncol=8)
  #browser()
  cvres <- cv_test3_c(xdec,y,lseq,xiv,bin0,1e-7,5)
  #cv0 <- cv_test1_c(x,y,lseq,bin0,1e-8,5)
  #cv2 <- cv_test2_c(x,y,lseq,bin0,1e-9,5)
  b0 <- cvres[[1]]
  b2 <- cvres[[4]]
  browser()
  res1 <- c(sum(abs(b0[1:k])>1e-10),sum(abs(b0[-(1:k)])>1e-10),sum((b0-b)^2),sum(abs(b2[1:k])>1e-10),sum(abs(b2[-(1:k)])>1e-10),sum((b2-b)^2))
  return(res1)
}

rr <- foreach(i=1:500,.packages = c("MASS","bar","mvtnorm"),.combine="rbind") %dopar% foo2(30,100,3,0.1)
rr <- foreach(i=1:500,.packages = c("MASS","bar","mvtnorm"),.combine="rbind") %dopar% foo3(30,100,3,0.1)

as.vector(itercd_c(x,y,p,n,0.2442053,bin0,1e-8)-itercd_c(x,y,p,n,0.2442053,ginv(t(x) %*% x+diag(rep(1e-5,p))) %*% t(x) %*% y,1e-8))

as.vector(barrep_c(x,y,0.2442053,bin0,1e-8))-as.vector(barrep_c(x,y,0.2442053,ginv(t(x) %*% x+diag(rep(100,p))) %*% t(x) %*% y,1e-8))

as.vector(barrep_c(x,y,0.2442053,ginv(t(x) %*% x+diag(rep(100,p))) %*% t(x) %*% y,1e-8))-as.vector(itercd_c(x,y,p,n,0.2442053,ginv(t(x) %*% x+diag(rep(100,p))) %*% t(x) %*% y,1e-8))

sapply(exp(seq(log(200),log(1e-5),length.out=100)),function(dd) max(as.vector(barrep_c(x,y,0.2442053,ginv(t(x) %*% x+diag(rep(dd,p))) %*% t(x) %*% y,1e-8))-as.vector(itercd_c(x,y,p,n,0.2442053,ginv(t(x) %*% x+diag(rep(dd,p))) %*% t(x) %*% y,1e-8))))


rr1 <- foreach(i=1:500,.packages = c("MASS","bar","mvtnorm"),.combine="rbind") %dopar% foo2(100,10,3,0.1)
rr2 <- foreach(i=1:100,.packages = c("MASS","bar","mvtnorm","RSpectra"),.combine="rbind") %dopar% foo3(100,10,3,0.1)
rr1 <- foreach(i=1:100,.packages = c("MASS","bar","mvtnorm"),.combine="rbind") %dopar% foo2(20,50,3,0.1)
rr2 <- foreach(i=1:100,.packages = c("MASS","bar","mvtnorm","RSpectra"),.combine="rbind") %dopar% foo3(20,50,3,0.1)

rr <- foreach(i=1:500,.packages = c("MASS","bar","mvtnorm"),.combine="rbind") %dopar% foo2(20,50,3,0.1)
colMeans(rr)

rr <- foreach(i=1:500,.packages = c("MASS","bar","mvtnorm"),.combine="rbind") %dopar% foo2(50,100,3,0.1)
colMeans(rr)

rr <- foreach(i=1:500,.packages = c("MASS","bar","mvtnorm"),.combine="rbind") %dopar% foo2(100,200,3,0.1)
colMeans(rr)
foot <- function(n,p){
  b <- c(rep(1,3),rep(0,p-3))
  x <- rmvnorm(n,rep(0,p),toeplitz(0.3^(1:p)))
  y <- x %*% b+t(rmvnorm(1,rep(0,n),diag(1,n)))*1
  A <- t(x) %*% x
  
  bin0 <- ginv(t(x) %*% x+diag(rep(1,p))) %*% t(x) %*% y
  bres <- summary(microbenchmark(itercd_c(x,y,dim(x)[2],dim(x)[1],0.1,bin0,1e-5),
                         barrep_c(x,y,0.1,bin0,1e-5),times=10L),unit="ms")
  bres$mean
}

tres <- mapply(foot,expand.grid(c(30,50,100),c(10,50,100,250,500,750,1000))[,1],expand.grid(c(30,50,100),c(10,50,100,250,500,750,1000))[,2])
tres <- cbind(expand.grid(c(30,50,100),c(10,50,100,250,500,750,1000)),t(tres))
colnames(tres) <- c("n","p","coordinate-wise update","ridge regression")
tdam <- melt(tres,id=c("n","p"))
colnames(tdam) <- c("n","p","Methods","time")
tdam$n <- as.character(tdam$n)
ggplot(tdam, aes(x=p, y=`time(milliseconds)`, colour=Methods, group=interaction(n,Methods),shape=n)) + 
  geom_point() + geom_line()

plog <- ggplot(tdam, aes(x=p, y=log(time), colour=Methods, group=interaction(n,Methods),shape=n)) + 
  geom_point() + geom_line()+scale_x_continuous("p") +
  scale_y_continuous("log(time) (milliseconds)")+theme(legend.position = "bottom")
pori <- ggplot(tdam, aes(x=p, y=time, colour=Methods, group=interaction(n,Methods),shape=n)) + 
  geom_point() + geom_line()+scale_x_continuous("p") +
  scale_y_continuous("time (milliseconds)")+theme(legend.position = "bottom")
  

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


barridge_c(A,B,n,p,1,0.1,1e-10)[[1]]-
barcd_c(A,B,n,p,1,0.1,1e-10)[[1]]

n  <- 1000
p <- 10
x=matrix(rnorm(n*p),n,p)
std <- .Call("standardize", x,F,PACKAGE="bar")
x <- std[[1]]
y=x %*% c(1,1,1,rep(0,p-3))+rnorm(n)
A <-t(x) %*% x
B <- t(x) %*% y

xp <- diag(sqrt(eigen(A)$values[1:10])) %*% t(eigen(A)$vector[,1:10])
max(t(xp) %*% xp -A)
barridge_c(A,B,n,p,1,0.1,1e-5)
barcd_c(A,B,n,p,1,0.1,1e-5)

microbenchmark(barridge_c(A,B,n,p,1,0.1,1e-5),barcd_c(A,B,n,p,1,0.1,1e-5))

system.time(barridge_c(A[1:3,1:3],B[1:3],n,3,1,0.1,1e-15))
microbenchmark(barcd_c(A,B,n,p,1,0.1,1e-15),times=50L)


microbenchmark(barcd_c(A,B,n,p,1,0.1,1e-5))

penreg_c(A,B,n,p,exp(seq(log(37.8),log(37.8*0.01),length.out=100)),rep(1,p),1,3.7,"scad",1e-4)

penreg_c(A,B,n,p,exp(seq(log(100),log(0.01),length.out=100)),rep(1,p),1,3.7,"scad",1e-4)

aa <- penreg_c(A,B,n,p,exp(seq(log(113),log(100*0.0001),length.out=100)),rep(1,p),1,3.7,"elnet",1e-5)
1-apply(aa[[1]],2,function(b) sum((y-x %*%b)^2))/ sum((y-mean(y))^2)

traint <- list(list(t(x[1:10,]) %*% x[1:10,],t(x[1:10,]) %*% y[1:10]),list(t(x[11:20,]) %*% x[11:20,],t(x[11:20,]) %*% y[11:20]))
testt <- list(list(t(x[11:20,]) %*% x[11:20,],t(x[11:20,]) %*% y[11:20]),list(t(x[1:10,]) %*% x[1:10,],t(x[1:10,]) %*% y[1:10]))

barcv_c(A, B, traint, testt, n, p, 1, c(10,0.9), 2, "cd", 1e-5)
barcv_c(A, B, traint, testt, n, p, c(1000,0.01), c(100,10,0.9,0.1), 2, "cd", 1e-5)


barcv_c(A, B, traint, testt, n, p, c(log(100),log(1),length.out=5), c(log(10.64156),log(10.64156*0.001),length.out=20), 2, "cd", 1e-5)
pencv_c(A, B, traint, testt, n, p, exp(seq(log(30.00574),log(30.00574*0.001),length.out=20)), rep(1,p), 2, 1, 3.7, "lasso", 1e-5)
pencv_c(A, B, traint, testt, n, p, exp(seq(log(30.00574),log(30.00574*0.01),length.out=20)), rep(1,p), 2, 1, 3.7, "scad", 1e-5)


penreg_c(A,B,n,p,5.253306,rep(1,p),1,3.7,"scad",1e-8)
at <- penreg(A,B,n,p,NULL,50,NULL,1,3.7,"scad")

plot(x=lv,y=opp[1,],xlim =c(max(lv),0),typ="l",ylim=c(-1.5,1.5))

for(i in 1:p){
  lines(x=lv,y=opp[i,])
}

st <-Sys.time()
1+1
et <- Sys.time()
et-st


load("C:/Users/commi/OneDrive - University of Missouri/PhD/Research/Recurrent Event Data with Broken Adaptive Ridge Regression/Rsim/bar/temptest.RData")
orit <- reggen(a[[2]])
cvcont <- cv.reggen(a[[2]],5)
btt <- pencv_c(orit$Omega,orit$Pv,cvcont$trainL,cvcont$testL,100,9,exp(seq(log(max(abs(orit$Pv))),log(max(abs(orit$Pv)*0.01)),length.out=200)),rep(1,9),5,0.1,3.7,"scad",1e-10)
btt$betaopt
pencv_c(orit$Omega,orit$Pv,cvcont$trainL,cvcont$testL,100,6,exp(seq(log(603),log(50),length.out=100)),rep(1,6),5,0.1,3.7,"lasso",1e-5)
penreg(orit$Omega,orit$Pv,100,6,NULL,100,NULL,1,3.7,"lasso",1e-8)



n  <- 100
p <- 10
x=matrix(rnorm(n*p,mean=0,sd=1),n,p)
y=x %*% c(1,1,1,rep(0,p-3))+rnorm(n)
A <-t(x) %*% x
B <- t(x) %*% y

As <- standard_c(A,B,1)[[1]]
Bs <- standard_c(A,B,1)[[2]]
std <- standard_c(A,B,1)[[3]]

Bs^2/diag(As)/4-((Bs-(As-diag(diag(As)))%*% solve(As+diag(0,p),Bs))^2/4/diag(As))

((B-(A-diag(diag(A)))%*% solve(A+diag(10000,p),B))^2/4/diag(A))

stdres <- barreg(A,B,n,p,xiv=exp(seq(log(1e-10),log(1e-10),length.out=10)),method="ridge")
nstdres <- barreg(A,B,n,p,xiv=exp(seq(log(100),log(0.01),length.out=10)),method="ridge",stand=F)

stdres <- barreg(A,B,n,p,method="ridge")
nstdres <- barreg(A,B,n,p,method="ridge",stand=F)
stdres[[1]][,1:50]-nstdres[[1]][,1:50]

b0 <- as.vector(solve(As+diag(1,p),Bs))

for(i in 1:50){
  b1 <- as.vector(barrdglupd_c(As, Bs, 35, b0))
  print(as.vector(b1),digits=7,width=6)
  b0 <- b1
}

foo <- function(x) x-as.vector(barrdglupd_c(As, Bs, 1, c(x,b0[-1])))[1]

reggen2 <- function(datac){
  datac <- datac[,order(unlist(datac["C",]))]
  Cv <- unlist(datac["C",])
  Zmat <- do.call(cbind,datac["Z",])
  olist <- unlist(datac[2,])
  tlist <- unlist(datac[3,])
  #browser()
  reggenint_c(Cv,Zmat,olist,tlist,length(tlist)*(rep(1:length(Cv),unlist(datac["m",]))-1)+1:length(tlist)-1)
  
  #Pv <- P_c(Zmat,tlist,olist,Cv)
  #return(1)
}
