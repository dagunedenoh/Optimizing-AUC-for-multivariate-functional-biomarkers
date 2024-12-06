#########################################################################################
# Methods
#########################################################################################

#Library load
library(pracma); library(mvtnorm) ; library(funData) ; library(MFPCA); library(pROC); library(fda.usc);library(glmnet)
library(Rsolnp) ; library(evmix); library(fda) ; library(statmod) ; library(pracma) ; library(splines2)
library(fdasrvf) ; library(fda) ; library(statmod) ; library(caret); library(bspline) ;library(splines)

#########################################################################################
# MFPCOC
#########################################################################################
ourM <- function(preDat, preDatD){
  
  #c <- rep(1,mplus)
  par(mfrow=c(1,4))
  r.1 <- c()
  r.2 <- c()
  preDat <- preDat[,-ind]
  preDatD <- preDatD[,-ind2]

  #for each fold
  for (i in 1:5){
    Y <- Y.new[-folds[[i]]]#trainY 222
    testY <- Y.new[folds[[i]]] #testY 56
    
    B.data <- preDat[-folds[[i]],] #train B.data from each fold
    D.data <- preDatD[-folds[[i]],] #train D.data
    test.B.data <- preDat[folds[[i]],] #test B.data
    test.D.data <- preDatD[folds[[i]],] #test D.data
    
    p1 <- funData(1:(ind-1),as.matrix(B.data)) 
    p2 <- funData(1:(ind2-1),as.matrix(D.data))
    
    #make multi Function Data
    d <- multiFunData(list(p1,p2))
    
    #MFPCA
    sp <- MFPCA(d, M=4,  uniExpansions = list(list(type = "uFPCA", 4),
                                              list(type = "uFPCA", 4)),fit = TRUE) 
    
    mplus <- which(c(summary(sp)[3,])>0.99)[1]-1
    sp <- MFPCA(d, M=mplus, uniExpansions = list(list(type = "uFPCA", 4),
                                                 list(type = "uFPCA", 4)),fit = TRUE) 

    #train fpc and reconstruct test fpc
    fpc1 <- (as.data.frame(sp$functions[[1]]@X)) #mplus*time
    fpc2 <- (as.data.frame(sp$functions[[2]]@X)) #mplus*time

    test_scores <- (as.matrix(test.B.data- colMeans(B.data))%*%t(as.matrix((fpc1)))) + (as.matrix(test.D.data- colMeans(D.data))%*%t(as.matrix((fpc2))))
    test_Y<- cbind(test_scores,(testY)) #test score matrix
    dat_Y<- cbind(sp$scores,(Y)) #train score matrix

    J<- function(c){
      
      Y_1 <- t(dat_Y[dat_Y[,dim(dat_Y)[2]]==1,-(mplus+1)]) 
      n_1 <- dim(Y_1)[2]
      Y_2 <- t(dat_Y[dat_Y[,dim(dat_Y)[2]]!=1,-(mplus+1)])
      n_2 <- dim(Y_2)[2]
      
      z_1 <- c%*%(Y_1) ; z_2 <- c%*%(Y_2) #composite scores
      
      #nonparametric approach based on AUC
      s=0
      for (n in 1:n_1){
        s2 = 0
        for (m in 1:n_2){
          s3 <- I(z_1[,n] > z_2[,m])
          s2 <- s2 + s3
        }
        s <- s + s2
      }
      
      J <- s/(n_1*n_2) 
      return(-J) #default is minimize. 
    }

    init <- ifelse(colMeans(sp$scores[which(Y==1),])>0,1,-1) #initial value
    res <- optim(par=c(init),fn=J) #optimize
    print(res$par)
    r.1[i] <- -res$value
    
    Y_1_ <- t(test_Y[test_Y[,dim(test_Y)[2]]==1,-(mplus+1)]) #diseased group = obstruct group / 3*74
    n_1_ <- dim(Y_1_)[2]
    Y_2_ <- t(test_Y[test_Y[,dim(test_Y)[2]]!=1,-(mplus+1)]) #healthy group / 3*204
    n_2_ <- dim(Y_2_)[2]
    
    #Composite score with optimal value
    z_1 <- res$par%*%(Y_1_) ; z_2 <- res$par%*%(Y_2_)

    s=0
    for (l in 1:n_1_){
      s2 = 0
      for (r in 1:n_2_){
        s3 <- I(z_1[,l] > z_2[,r])
        s2 <- s2 + s3
      }
      s <- s + s2
    }
    v <- s/(n_1_*n_2_) 
    r.2[i] <- v
  }
  return(list(r.1, r.2))
}

#########################################################################################
# SPFCOC in multivariate setting
#########################################################################################
ourM2 <- function(preDat, preDatD){
  
  B.data.new <- preDat[,-ind]
  D.data.new <- preDatD[,-ind2]
  r.1 <- c()
  r.2 <- c()
  
  for (i in 1:5){
    
    Y <- Y.new[-folds[[i]]]#trainY
    testY <- Y.new[folds[[i]]] #testY
    
    B.data <- B.data.new[-folds[[i]],] #train B.data
    D.data <- D.data.new[-folds[[i]],] #train D.data
    test.B.data <- B.data.new[folds[[i]],] #test B.data
    test.D.data <- D.data.new[folds[[i]],] #test D.data
    
    #Bspline
    Bspline2 <- create.bspline.basis(rangeval = c(1,(ind-1)),nbasis=8)
    Bspline <- create.bspline.basis(rangeval = c(1,(ind2-1)),nbasis=8)
    
    #To functional data
    B.fd <- Data2fd(argvals = c(1:(ind-1)), y=t(B.data),basisobj = Bspline2)
    D.fd <- Data2fd(argvals = c(1:(ind2-1)), y=t(D.data),basisobj = Bspline)
    test.B.fd <- Data2fd(argvals = c(1:(ind-1)), y=t(test.B.data),basisobj = Bspline2)
    test.D.fd <- Data2fd(argvals = c(1:(ind2-1)), y=t(test.D.data),basisobj = Bspline)
    
    #FPC
    a <- pca.fd(B.fd, nharm=8) #PC
    a_n <- (which( cumsum(a$varprop) > 0.99))[1]-1
    a <- pca.fd(B.fd, nharm=a_n) #PC
    
    b <- pca.fd(D.fd, nharm=8) #PC
    b_n <- (which(cumsum(b$varprop) > 0.99))[1]-1
    b <- pca.fd(D.fd, nharm=b_n) #PC
    
    #Stack FPC score matrix
    dat_Y<- cbind(a$scores, b$scores,(Y))
    mat1 <-  bs(x = 1:59,  df=8,  intercept = FALSE)
    mat2 <-  bs(x = 1:40,  df=8,  intercept = FALSE)
    
    #Reconstruct score matrix
    test1 <- t(mat1%*%(as.matrix(test.B.fd$coefs) - as.vector(a$meanfd$coefs)))%*%mat1%*%(as.matrix((a$harmonics$coefs)))
    test2 <- t(mat2%*%(as.matrix(test.D.fd$coefs) - as.vector(b$meanfd$coefs)))%*%mat2%*%(as.matrix((b$harmonics$coefs)))
    
    #Test FPC score matrix
    test_Y <- (cbind(test1, test2, testY))
    
    J<- function(c){
      
      Y_1 <- t(dat_Y[dat_Y[,dim(dat_Y)[2]]==1,-(a_n+b_n+1)]) #diseased group
      n_1 <- dim(Y_1)[2]
      Y_2 <- t(dat_Y[dat_Y[,dim(dat_Y)[2]]!=1,-(a_n+b_n+1)]) #healthy group
      n_2 <- dim(Y_2)[2]
      
      z_1 <- c%*%(Y_1) ; z_2 <- c%*%(Y_2) #composite score
      
      s=0
      for (n in 1:n_1){
        s2 = 0
        for (m in 1:n_2){
          s3 <- I(z_1[,n] > z_2[,m])
          s2 <- s2 + s3
        }
        s <- s + s2
      }
      
      J <- s/(n_1*n_2)
      return(-J)
    }
    init <- ifelse(colMeans(cbind(a$scores,b$scores)[which(Y==1),])>0,1,-1)
    res <- optim(par=c(init),fn=J) #optimize
    r.1[i] <- -res$value #optimal coefficient
    
    Y_1_ <- t(test_Y[test_Y[,dim(test_Y)[2]]==1,-(a_n+b_n+1)]) #diseased group
    n_1_ <- dim(Y_1_)[2]
    Y_2_ <- t(test_Y[test_Y[,dim(test_Y)[2]]!=1,-(a_n+b_n+1)]) #healthy group
    n_2_ <- dim(Y_2_)[2]
    
    z_1 <- res$par%*%(Y_1_) ; z_2 <- res$par%*%(Y_2_) #Copmosite score with optimal coefficient

    s=0
    for (l in 1:n_1_){
      s2 = 0
      for (r in 1:n_2_){
        s3 <- I(z_1[,l] > z_2[,r])
        s2 <- s2 + s3
      }
      s <- s + s2
    }
    v <- s/(n_1_*n_2_)
    r.2[i] <- v
  }
  return(list(r.1, r.2))
}

#########################################################################################
# SFPCOC in univariate setting
#########################################################################################
ourM3 <- function(preDat){
  
  B.data.new <- preDat[,-ind]
  r.1 <- c()
  r.2 <- c()
  
  for (i in 1:5){
    print(i)
    
    Y <- Y.new[-folds[[i]]]#trainY
    testY <- Y.new[folds[[i]]] #testY
    
    B.data <- B.data.new[-folds[[i]],] #train B.data
    test.B.data <- B.data.new[folds[[i]],] #test B.data

    #Bspline
    Bspline2 <- create.bspline.basis(rangeval = c(1,(ind-1)),nbasis=8)
    
    #Functional Data : train
    B.fd <- Data2fd(argvals = c(1:(ind-1)), y=t(B.data),basisobj = Bspline2)
    test.B.fd <-Data2fd(argvals = c(1:(ind-1)), y=t(test.B.data),basisobj = Bspline2)
    
    a <- pca.fd(B.fd, nharm=8) #PC
    a_n <- (which( cumsum(a$varprop) > 0.99))[1]-1
    a <- pca.fd(B.fd, nharm=a_n) #PC
    
    mat1 <-  bs(x = 1:59,  df=8,  intercept = FALSE)
    mat2 <-  bs(x = 1:40,  df=8,  intercept = FALSE)
    
    #calculate test FPC score
    test_scores <- t(mat1%*%(as.matrix(test.B.fd$coefs) - as.vector(a$meanfd$coefs)))%*%mat1%*%(as.matrix((a$harmonics$coefs)))
    
    dat_Y<- cbind(a$scores,(Y))
    test_Y <- cbind(test_scores, testY)
    
    J<- function(c){
      
      Y_1 <- t(dat_Y[dat_Y[,dim(dat_Y)[2]]==1,-(a_n+1)]) #diseased group = obstruct group / 3*74
      n_1 <- dim(Y_1)[2]
      Y_2 <- t(dat_Y[dat_Y[,dim(dat_Y)[2]]!=1,-(a_n+1)]) #healthy group / 3*204
      n_2 <- dim(Y_2)[2]
      
      z_1 <- c%*%(Y_1) ; z_2 <- c%*%(Y_2) #composite score
      
      s=0
      for (n in 1:n_1){
        s2 = 0
        for (m in 1:n_2){
          s3 <- I(z_1[,n] > z_2[,m])
          s2 <- s2 + s3
        }
        s <- s + s2
      }
      J <- s/(n_1*n_2)
      return(-J) #default is minimize.
    }
    res <- optim(par=c(rep(1,a_n)),fn=J) #optimize
    r.1[i] <- -res$value #optimal coefficient
    
    #Test data
    Y_1_ <- t(test_Y[test_Y[,dim(test_Y)[2]]==1,-(a_n+1)]) #diseased group = obstruct group / 3*74
    n_1_ <- dim(Y_1_)[2]
    Y_2_ <- t(test_Y[test_Y[,dim(test_Y)[2]]!=1,-(a_n+1)]) #healthy group / 3*204
    n_2_ <- dim(Y_2_)[2]
    
    #Composite score with optimal coefficient
    z_1 <- res$par%*%(Y_1_) ; z_2 <- res$par%*%(Y_2_)
    
    #maximize Mann- Whitney statistics
    s=0
    for (l in 1:n_1_){
      s2 = 0
      for (r in 1:n_2_){
        s3 <- I(z_1[,l] > z_2[,r])
        s2 <- s2 + s3
      }
      s <- s + s2
    }
    v <- s/(n_1_*n_2_)
    r.2[i] <- v
  }
  return(list(r.1, r.2))
}

#########################################################################################
# FPCLR (logit FD) in multivariate setting
#########################################################################################

logitFDF2 <- function(B.data.new, D.data.new, id){
  
  library(logitFD)
  B.data.new <- preDat[,-ind]
  D.data.new <- preDatD[,-ind2]
  true_Y <- preDat[,ind]
  true_Y <- ifelse(true_Y==-1,0,1)
  r.1 <- c()
  r.2 <- c()
  for(i in 1:5){
    print(i)
    
    Y.new <- ifelse(Y.new==1,1,0)
    Y <- Y.new[-folds[[i]]]#trainY
    testY <- Y.new[folds[[i]]] #testY
    
    B.data <- B.data.new[-folds[[i]],] #train B.data
    D.data <- D.data.new[-folds[[i]],] #train B.data
    test.B.data <- B.data.new[folds[[i]],] #test B.data
    test.D.data <- D.data.new[folds[[i]],] #test B.data
    
    #Bspline
    Bspline2 <- create.bspline.basis(rangeval = c(1,(ind-1)),nbasis=8)
    Bspline <- create.bspline.basis(rangeval = c(1,(ind2-1)),nbasis=8)
    
    #Functional Data : train
    B.fd <- Data2fd(argvals = c(1:(ind-1)), y=t(B.data),basisobj = Bspline2)
    D.fd <- Data2fd(argvals = c(1:(ind2-1)), y=t(D.data),basisobj = Bspline)
    
    res <- logitFD.pc(Response=Y, FDobj = list(B.fd,D.fd), ncomp = c(8,8), nonFDvars = NULL)
    a_n <- which( res$PC.variance[[1]]$`% Cum.Prop.Var`>99)[1]-1
    b_n <- which( res$PC.variance[[2]]$`% Cum.Prop.Var`>99)[1]-1
    res <- logitFD.pc(Response=Y, FDobj = list(B.fd,D.fd), ncomp = c(a_n,b_n), nonFDvars = NULL)
    
    
    a <- pca.fd(B.fd, nharm=a_n) #PC for test data
    b <- pca.fd(D.fd, nharm=b_n) #PC for train data
    
    r.1[i] <- res$ROC$auc/100 #train
    
    #Functional data : test
    test.B.fd <- Data2fd(argvals = c(1:(ind-1)), y=t(test.B.data),basisobj = Bspline2)
    test.D.fd <- Data2fd(argvals = c(1:(ind2-1)), y=t(test.D.data),basisobj = Bspline)
    
    mat1 <-  bs(x = 1:59,  df=8,  intercept = FALSE)
    mat2 <-  bs(x = 1:40,  df=8,  intercept = FALSE)
    
    a_scores <- t(mat1%*%(as.matrix(test.B.fd$coefs) - as.vector(a$meanfd$coefs)))%*%mat1%*%(as.matrix((a$harmonics$coefs)))
    b_scores <- t(mat2%*%(as.matrix(test.D.fd$coefs) - as.vector(b$meanfd$coefs)))%*%mat2%*%(as.matrix((b$harmonics$coefs)))
    
    #Test processs
    scoreMat <- cbind(1,a_scores, b_scores)
    ks <- (scoreMat%*%coef(res$glm.fit))
    res2 <- as.vector(ifelse((exp(ks)/(exp(ks)+1))>0.5,1,0))
    res2 <- replace(res2,is.na(res2),1)
    r.2[i]<- auc(roc(testY,res2)) #test
  }
  return(list(r.1,r.2))
}

#########################################################################################
# FPCLR (logit FD) in univariate setting
#########################################################################################

logitFDF <- function(B.data.new, id){
  
  B.data.new <- preDat[,-ind]
  library(logitFD)
  true_Y <- preDat[,ind]
  true_Y <- ifelse(true_Y==-1,0,1)
  r.1 <- c()
  r.2 <- c()
  
  for(i in 1:5){
    Y.new <- ifelse(Y.new==1,1,0)
    Y <- Y.new[-folds[[i]]]#trainY
    testY <- Y.new[folds[[i]]] #testY
    
    B.data <- B.data.new[-folds[[i]],] #train B.data
    test.B.data <- B.data.new[folds[[i]],] #test B.data
    #Bspline
    Bspline2 <- create.bspline.basis(rangeval = c(1,(ind-1)),nbasis=8)
    
    #Functional Data : train
    B.fd <- Data2fd(argvals = c(1:(ind-1)), y=t(B.data),basisobj = Bspline2)
    
    res <- logitFD.pc(Response=Y, FDobj = list(B.fd), ncomp = c(8), nonFDvars = NULL)
    a_n <- which( res$PC.variance[[1]]$`% Cum.Prop.Var`>99)[1]-1
    res <- logitFD.pc(Response=Y, FDobj = list(B.fd), ncomp = c(a_n), nonFDvars = NULL)
    
    r.1[i] <- res$ROC$auc/100 #train
    
    #Functional data : test
    a <- pca.fd(B.fd, nharm=a_n) #train data
    test.B.fd <- Data2fd(argvals = c(1:(ind-1)), y=t(test.B.data),basisobj = Bspline2)
    
    mat1 <-  bs(x = 1:59,  df=8,  intercept = FALSE)
    a_scores <- t(mat1%*%(as.matrix(test.B.fd$coefs) - as.vector(a$meanfd$coefs)))%*%mat1%*%(as.matrix((a$harmonics$coefs)))
    
    #Test processs
    scoreMat <- cbind(1,a_scores)
    ks <- (scoreMat%*%coef(res$glm.fit))
    res2 <- as.vector(ifelse((exp(ks)/(exp(ks)+1))>0.5,1,0))
    res2 <- replace(res2,is.na(res2),1)
    r.2[i]<- auc(roc(testY,res2)) #test
  }
  return(list(r.1,r.2))
}

########################################################################################
# Functional Classification (Logistic) in multivariate setting
#########################################################################################

fLogistic2 <- function(B.data.new, D.data.new, Y.new, id){
  r.1 <- c()
  r.2 <- c()
  
  for(i in 1:5){
    
    Y <- Y.new[-folds[[i]]]#trainY
    testY <- Y.new[folds[[i]]] #testY
    
    B.data <- B.data.new[-folds[[i]],] #train B.data
    D.data <- D.data.new[-folds[[i]],] #train B.data
    test.B.data <- B.data.new[folds[[i]],] #test B.data
    test.D.data <- D.data.new[folds[[i]],] #test B.data
    
    x <-fdata(B.data); x2 <-fdata(D.data)
    
    ldat <- ldata("df"=data.frame(y=factor(Y)),"x"=x, "x2"=x2)
    a1 <- classif.glm(y ~ x + x2, data = ldat, prob=0.5)
    
    p1 <- as.numeric(predict(a1,ldat))#predict
    p1 <- ifelse(p1==2,1,0)
    p2 <- ifelse(ldat$df$y==1,1,0)
    ROC <- roc(p2,p1)
    r.1[i] <-auc(ROC)
    
    x <-fdata(test.B.data)
    x2 <-fdata(test.D.data)
    
    newldat <- ldata("df"=data.frame(y=testY),"x"=x, "x2"=x2)
    p1 <- as.numeric(predict(a1,newldat))#predict
    p1 <- ifelse(p1==2,1,0)
    p2 <- ifelse(newldat$df$y==1,1,0)
    ROC <- roc(p2,p1)
    #print(auc(ROC))
    r.2[i] <-auc(ROC)
  }
  return(list(r.1,r.2))
}

#########################################################################################
# Functional Classification in univariate setting
#########################################################################################

fLogistic <- function(B.data.new, Y.new, id){
  r.1 <- c()
  r.2 <- c()
  for(i in 1:5){    
    Y <- Y.new[-folds[[i]]]#trainY
    testY <- Y.new[folds[[i]]] #testY
    
    B.data <- B.data.new[-folds[[i]],] #train B.data
    test.B.data <- B.data.new[folds[[i]],] #test B.data
    x <-fdata(B.data)
    
    ldat <- ldata("df"=data.frame(y=factor(Y)),"x"=x)
    a1 <- classif.glm(y ~ x, data = ldat, prob=0.5)
    
    p1 <- as.numeric(predict(a1,ldat))#predict
    p1 <- ifelse(p1==2,1,0)
    p2 <- ifelse(ldat$df$y==1,1,0)
    ROC <- roc(p2,p1)
    r.1[i] <-auc(ROC)
    
    x <-fdata(test.B.data)
    
    newldat <- ldata("df"=data.frame(y=testY),"x"=x)
    p1 <- as.numeric(predict(a1,newldat))#predict
    p1 <- ifelse(p1==2,1,0)
    p2 <- ifelse(newldat$df$y==1,1,0)
    ROC <- roc(p2,p1)
    r.2[i] <-auc(ROC)
  }
  return(list(r.1,r.2))
}
