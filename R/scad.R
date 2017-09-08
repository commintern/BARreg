penreg <- function(A,B,n,p,lambdav=NULL,nlambda=50,weight=NULL, alpha=1,gamma=3.7,method=c("scad","lasso","alasso"),abstol=1e-8,ifc=T){
  method <- match.arg(method)
  if(!ifc){
    B <- t(A) %*% B
    A <- t(A) %*% A
  }
  if(is.null(lambdav)){
    lambda_max <- max(abs(B))
    lambda.min=ifelse(n>p,.001,.05)
    lambdav <- exp(seq(log(lambda_max),log(lambda_max*lambda.min),length.out = nlambda+1))[-1]
  }
  if(method=="lasso"){
    res <- penreg_c(A,B,n,p,lambdav,rep(1,p),alpha=1,gamma=3.7,method,abstol)
  } 
  if(method=="alasso"){
    weight <- 1/abs(solve(A+diag(nrow=p),B))
    res <- penreg_c(A,B,n,p,lambdav,weight,alpha=1,gamma=3.7,method,abstol)
  }
  if(method=="scad"){
    res <- penreg_c(A,B,n,p,lambdav,rep(1,p),alpha=1,gamma=3.7,method,abstol)
  }
  res
}

cv.genpen <- function(oricon, cvcon, n, p, nfold = 5, nlambda = 50, abstol = 1e-08, method = c("lasso", "alasso", "scad"), alsp = 0.01,lambda.min.ratio=1e-3,ifc = T) {
  method <- match.arg(method)
  #oricon <- reggen(datac)
  lambda_max <- max(abs(oricon$Pv))
  lambdav <- exp(seq(log(lambda_max*(1)), log(lambda_max * lambda.min.ratio), length.out = nlambda))
  #cvcon <- cv.reggen(datac, nfold)
  
  if (method == "lasso") {
    strTime <- Sys.time()
    res <- pencv_c(oricon$Omega, oricon$Pv, cvcon$trainL, cvcon$testL, n, p, lambdav, rep(1, p), nfold, alpha = 1, gamma = 3.7, method, abstol)
    endTime <- Sys.time()
  }
  if (method == "alasso") {
    weight <- as.vector(1/abs(solve(oricon$Omega + alsp*diag(nrow = p), oricon$Pv)))
    strTime <- Sys.time()
    res <- pencv_c(oricon$Omega, oricon$Pv, cvcon$trainL, cvcon$testL, n, p, lambdav, weight, nfold, alpha = 1, gamma = 3.7, method, abstol)
    endTime <- Sys.time()
  }
  if (method == "scad") {
    strTime <- Sys.time()
    res <- pencv_c(oricon$Omega, oricon$Pv, cvcon$trainL, cvcon$testL, n, p, lambdav, rep(1, p), nfold, alpha = 1, gamma = 3.7, method, abstol)
    endTime <- Sys.time()
  }
  if(length(which(res[[3]][[3]]==10000)>0)) warning(paste(method," fails to converge"))
  c(res,endTime-strTime)
}