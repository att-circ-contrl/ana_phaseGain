#! /usr/bin/Rscript

library("jsonlite")
library("glmnet")
library("R.matlab")

# testing

# unfurl the command line argument
args = commandArgs(TRUE)
tmpname = args[1]
inpath = paste(tmpname,"_in.mat",sep="")
outpath = paste(tmpname,"_out.mat",sep="")

dat <- readMat(inpath)

# extract the data
family = dat[["family"]]
X = dat[["X"]]
y = dat[["y"]]
nrand = dat[["nrand"]]
nfold = dat[["nfold"]]

if (!is.null(dat[["useMinLambda"]])){
  useMinLambda = dat[["useMinLambda"]]
}else{
  useMinLambda = 0
}

# ignore rows with missing values
good = rowSums(is.na(X))==0;
X = X[good,]
y = y[good]

sz = dim(X)

# bootstraped CV
Bmin = matrix(0,nrand,sz[2]+1)
B1se = matrix(0,nrand,sz[2]+1)
for (ir in 1:nrand){
  print(ir)
  idx = sample.int(sz[1], size = sz[1], replace = TRUE)
  X2 = X[idx,]
  y2 = y[idx]

  fit = cv.glmnet(X2,y2,family=family,alpha=0.95,nfolds=nfold,keep=TRUE) #a=0.95, runyan 2017 nature

  B1se[ir,] = matrix(coef(fit,fit$lambda.1se))
  Bmin[ir,] = matrix(coef(fit,fit$lambda.min))
}

# prep output and save
if (useMinLambda){
  Ball = Bmin
}else{
  Ball = B1se
}

B = colMeans(Ball)

dummy.call <- t(c("x", "y", "family", "options"))
writeMat(outpath,Ball=Ball,Bmin=Bmin,B1se=B1se,B=B)
