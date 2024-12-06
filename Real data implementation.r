#######################################################################################
# * MFPCOC for renal studies
#######################################################################################

library(Rsolnp) ; library(evmix) ; library(tidyverse); library(splitTools)
library(fda) ; library(statmod) ; library(pracma) ; library(splines2)
library(fdasrvf) ; library(fda) ; library(statmod) ; library(caret); library(bspline) ;library(splines)

source("0827 Method_Real.R")
set.seed(613)

# 1. Data load 
data <- subset(data, binary_consensus!=2)
data$binary_consensus[data$binary_consensus==0] <- -1 
patient <- data[1:2, 1] #age, gender, AT, ED, RH, binary_consensus, Kidney

n <- 278
B <- data[1:n,11:69]
D <- data[1:n,70:dim(data)[2]]
Y <- data$binary_consensus #278*1 scale should be changed {-1,1}

B.dat <- cbind (patient, B)
B.data <- (B.dat[1:n,-1])
B.data <- as.matrix(B.data) #time point * patient number #278*59

D.dat <- cbind(patient,D)
D.data <- (D.dat[1:n,-1])
D.data <- as.matrix(D.data)

numbersB <- round(as.integer(gsub("B", "", colnames(B.data)))/60)
numbersD <- round(as.integer(gsub("D", "", colnames(D.data)))/60)

# 2. Calculate integrated variabilty

test <- matrix(c(1,2,3,4,5,2,4,6,8,10), nrow=2, ncol=5)
test_m <- colMeans(test)

ss=0
for (i in 1:dim(B.data)[1]){
  s <- trapz(1:59,(t(B.data[i,]) - colMeans(B.data))^2)
  ss <- s + ss
}

ss2=0
for (i in 1:dim(D.data)[1]){
  s <- trapz(1:40,(t(D.data[i,]) - colMeans(D.data))^2)
  ss2 <- s + ss2
}

c <- ss/ss2

C.data <- as.data.frame(cbind(B.data/1000, Y))

delta <- which(C.data$Y==1)
delta2 <- which(C.data$Y==-1)

Y.new<- C.data[,dim(C.data)[2]]
B.data.new <- C.data[,-dim(C.data)[2]]

Obs <- (B.data.new[delta,]) #time * number of curves : 59*74
Normal <- (B.data.new[delta2,]) # time * number of curves : 59*204
B.data.new_ <- cbind(rbind(Obs, Normal), c(rep(1,length(delta)), rep(-1,length(delta2))))
B.data.new_ <- as.matrix(B.data.new_)
matplot(t(B.data.new_[,-60]), type="l",main="1st renogram curve") #Plot

# D data scale
B.dat <- cbind (patient, B)
B.data <- (B.dat[1:n,-1])
B.data <- as.matrix(B.data) #time point * patient number #278*59

D.dat <- cbind(patient,D)
D.data <- (D.dat[1:n,-1])
D.data <- as.matrix(D.data)/1000
D.data <- D.data*sqrt(c)

C.data <- as.data.frame(cbind(D.data, Y))
delta <- which(C.data$Y==1)
delta2 <- which(C.data$Y==-1)
Y.new<- C.data[,dim(C.data)[2]]
D.data.new <- C.data[,-dim(C.data)[2]]
Obs <- (D.data.new[delta,]) #time * number of curves : 59*74
Normal <- (D.data.new[delta2,]) # time * number of curves : 59*204
D.data.new_ <- cbind(rbind(Obs, Normal), c(rep(1,length(delta)), rep(-1,length(delta2))))
D.data.new_ <- as.matrix(D.data.new_)

#final data
preDat <- B.data.new_
preDatD <- D.data.new_

#Basis
x = seq(1,dim(preDat)[2]-1)
ind <- dim(preDat)[2] #60
ind2 <- dim(preDatD)[2] #41
true_Y = t(preDat[,ind]) #same
Y1 = t(preDat[,-ind]) #59*278
Y2 = t(preDatD[,-ind2]) #40*278
Y <- rbind(Y1,Y2) #99*278
colMeans(B.data.new_[c(which(true_Y==1)),-60])

basisObj <- function(basis){
  if(basis == "FPC"){
    mat <- create.pc.basis(fdata(t(Y)), 1:3)
    plot(mat$basis, col=3)
    mat <- t(mat$basis$data) #59*3
  }
  if(basis=="BsplineTrue"){
    mat <-  bs(x = x,  df=10,  intercept = TRUE) #basis 개수 * time point 개수
    matplot(t(mat),type='l')
  }
  if(basis=="BsplineFalse"){
    mat <-  bs(x = x,  df=10,  intercept = FALSE) #basis 개수 * time point 개수
    matplot(t(mat),type='l')
  }
  return(mat)
}


################################################
# Implementation
################################################

preDat <- B.data.new_
preDatD <- D.data.new_

row.names(preDat) <- row.names(preDatD)
x <- as.double(rownames(preDat))

ind = dim(preDat)[2]
ind2 = dim(preDatD)[2]
Y.new <- true_Y <- preDat[,ind]

colnames(preDat)[ind] <- "Y" ; colnames(preDatD)[ind2] <- "Y"

#5-fold which is used in function M1, M2, M3, Flog, Flog2, LFd, LFd2
dis <- which(Y.new==1)
nondis <- which(Y.new==-1)
fold1 <- c(dis[1:15], nondis[1:41])
fold2 <- c(dis[16:30], nondis[42:82])
fold3 <- c(dis[31:45], nondis[83:123])
fold4 <- c(dis[46:60], nondis[124:164])
fold5 <- c(dis[61:74], nondis[165:204])
folds <- list(fold1, fold2, fold3, fold4, fold5)

################################################
# 1. ourM = MFPCOC
# 2. ourM2 = SFPCOC in multivariate setting
# 3. ourM3 = SFPCOC in univariate setting
################################################

M1 <- ourM(preDat[,-ind], preDatD[,-ind2])
M2 <- ourM2(preDat, preDatD)
M3 <- ourM3(preDat)

total[iter] <- mean(M1[[1]])
total2[iter] <- mean(M1[[2]])
totalB[iter] <- mean(M2[[1]])
total2B[iter] <- mean(M2[[2]])
totalC[iter] <- mean(M3[[1]])
total2C[iter] <- mean(M3[[2]])

############################
# F logistic : Separate
############################
Y.new <- true_Y
Flog <- fLogistic(preDat[,-ind],Y.new,id) #univariate setting
Flog2 <- fLogistic2(preDat[,-ind], preDatD[,-ind2],Y.new,id) #multivariate setting

resFlog3[iter] <- mean(Flog[[1]])
resFlog4[iter] <- mean(Flog[[2]])
resFlog5[iter] <- mean(Flog2[[1]])
resFlog6[iter] <- mean(Flog2[[2]])

############################
# logit FD 
############################

LFd <- logitFDF(preDat[,-ind],id)  #univariate setting
LFd2 <- logitFDF2(preDat[,-ind],preDatD[,-ind2],id) #multivariate setting
resFD1[iter] <- mean(LFd[[1]])
resFD2[iter] <- mean(LFd[[2]])
resFD3[iter] <- mean(LFd2[[1]])
resFD4[iter] <- mean(LFd2[[2]])

print(paste(round(mean(total2C),3), round(mean(resFlog4),3), round(mean(resFD2),3)))
print(paste(round(mean(total2),3), round(mean(total2B),3), round(mean(resFlog6),3), round(mean(resFD4),3)))
  

#Mean AUC
print(paste(round(mean(M1[[2]]),3), round(mean(M2[[2]]),3), round(mean(M3[[2]]),3)))
print(paste(round(mean(Flog[[2]]),3), round(mean(Flog2[[2]]),3), round(mean(LFd[[2]]),3), round(mean(LFd2[[2]]),3)))

#Standard deviation between 5-fold CV
print(paste(round(var(M1[[2]]),3), round(var(M2[[1]]),3), round(var(M3[[1]]),3)))
print(paste(round(var(Flog[[2]]),3), round(var(Flog2[[2]]),3), round(var(LFd[[2]]),3), round(var(LFd2[[2]]),3)))