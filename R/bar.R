# barreg <- function(X,Y,xki=1,lambdav=NULL,nlambda=100,abstol=1e-10,ifC=T，method=c("ridge")){
#   if (!requireNamespace("ncvreg", quietly = TRUE)) {
#     stop("ncvreg needed for this function to work. Please install it.",
#          call. = FALSE)
#   }
#   iter <- 0
#   std <- .Call("standardize", X,PACKAGE="ncvreg")
#   XX <- std[[1]]
#   center <- std[[2]]
#   scale <- std[[3]]
#   XX <- t(t(XX)+center/scale)
#   nz <- which(scale > 1e-6)
#   if (length(nz) != ncol(XX)) XX <- XX[ ,nz, drop=FALSE]
#   YY <- Y
#   betao<-solve(t(XX)%*%XX+xki*diag(nrow=length(XX[1,])))%*%t(XX)%*%YY
#   beta0 <- betao
#   if(is.null(lambdav)){
#     lambda_max <- max((t(XX) %*% Y)^2 / diag(t(XX) %*% XX) /4 )
#     lambdav <- exp(seq(log(lambda_max*0.001),log(lambda_max),length.out = nlambda))[-1]
#   }
#   lambdav <- sort(lambdav,decreasing=T)
#   if(ifC){
#     beta1 = barridge_c(t(XX) %*% XX,t(XX) %*%YY,dim(XX)[1],dim(XX)[2],exp(seq(log(100),log(0.01),length.out = 5)),lambdav,abstol)
#   } else {
#     repeat{
#       beta1 <- updatefun(XX,YY,lambdav[1],beta0)
#       if(mean(abs((beta1-beta0))) < abstol) break else beta0 <- beta1
#       iter <- iter +1 
#       if(iter>10000) {
#         cat("Fail to converge")
#         break
#       }
#     }
#   }
# 
#   beta1 <- beta1/scale[nz]
#   beta1 <- ifelse(beta1<1e-10,0,beta1)
#   list(beta=beta1,xki=xki,lambda=lambdav)
# }

barreg <- function(A,B,n,p,xiv=NULL,lambdav=NULL,nlambda=50,abstol=1e-8,method=c("ridge","cd"),lambda.min.ratio=0.001,ifc = T,stand=T){
  method <- match.arg(method)
  if(!ifc){
    B <- t(A) %*% B
    A <- t(A) %*% A
  }
  if(stand){
    temp <- standard_c(A,B,1)
    A <- temp[[1]]
    B <- temp[[2]]
    std <- as.vector(temp[[3]])
  }
  
  if(is.null(lambdav)){
    lambda_max <- max(B^2 / diag(A) /4 )
    lambdav <- exp(seq(log(lambda_max),log(lambda_max*lambda.min.ratio),length.out = nlambda))
  }
  lambdav <- sort(lambdav,decreasing=T)
  if(is.null(xiv)){
    xiv <- exp(seq(log(100), log(1e-4), length.out = 10))
  }
  if(method=="ridge"){
    res <- barridge_c(A,B,n,p,xiv,lambdav,abstol)
  } else {
    res <- barcd_c(A,B,n,p,xiv,lambdav,abstol)
  }
  if(stand){
    res[[1]] <- res[[1]] * std
  }
  
  res
}


cv.barreg <- function(oricon, cvcon, n, p, xiv=NULL,lambdav=NULL,nfold = 5, nlambda = 50, abstol = 1e-08, method = c("ridge", "cd"), lambda.min.ratio=ifelse(n>p,1e-3,0.05),stand=T,ifc = T) {
  method <- match.arg(method)
  #oricon <- reggen(datac)
  if(is.null(lambdav)){
    lambda_max <- max(oricon$Pv^2/diag(oricon$Omega)/4)
    lambdav <- exp(seq(log(lambda_max), log(lambda_max * lambda.min.ratio), length.out = nlambda))
  }
  if(is.null(xiv)){
    xiv <- exp(seq(log(10), log(1e-4), length.out = 5))
  }
  
  #cvcon <- cv.reggen(datac, nfold)
  strTime <- Sys.time()
  res <- barcv_c(oricon$Omega, oricon$Pv, cvcon$trainL, cvcon$testL, n, p, xiv, lambdav, nfold, method, abstol,0)
  endTime <- Sys.time()
  if(length(which(res[[3]][[3]]==10000)>0)) warning(paste(method," fails to converge"))
  c(res,endTime-strTime)
}
