reggen<- function(datac){
  datac <- datac[,order(unlist(datac["C",]))]
  Cv <- unlist(datac["C",])
  Zmat <- do.call(cbind,datac["Z",])
  olist <- unlist(datac[2,])
  tlist <- unlist(datac[3,])
  #browser()
  res <- reggenint_c(Cv,Zmat,olist,tlist,length(tlist)*(rep(1:length(Cv),unlist(datac["m",]))-1)+1:length(tlist)-1)
  list(Omega=res[[1]]/length(Cv),Pv=res[[2]]/length(Cv))
}


cv.reggen <- function(datac,nfold){
  no <- length(datac[1,])
  splitno <- split(1:no,1:no %% nfold)
  testL <- lapply(1:nfold,function(i) reggen(datac[,splitno[[i]]]))
  trainL <- lapply(1:nfold,function(i) reggen(datac[,-splitno[[i]]]))
  list(testL=testL,trainL=trainL)
}