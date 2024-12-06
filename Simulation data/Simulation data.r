#CORRELATED DATA
library(pracma);library(mvtnorm);library(funData); library(MFPCA)
library(Rsolnp) ; library(evmix); library(fda) ; library(statmod) ; library(pracma) ; library(splines2)
library(fdasrvf) ; library(fda) ; library(statmod) ; library(caret); library(bspline) ;library(splines)
library(pROC); library(fda.usc);library(glmnet)

########################################################################
# Simulation Data
# num : number of subjects
# time : time domain for 1st funcitonal object
# time2 : time domain for 2nd funcitonal object
# iter : allocate separate seed for different data among 300 data
########################################################################

random_numbers_with_sum <- function(n,iter) {
  n1 <- n*0.75
  n2 <- n*0.25
  
  return(c(n1,n2))
}

dataC <- function(num, time, time2,iter){
  set.seed(iter)
  
  num_sort = sort(random_numbers_with_sum(num,iter))
  print(num_sort)
  n1 = round(num_sort[1]*0.5)
  n2 = round(num_sort[1]*0.5)
  n3 = round(num_sort[2]*0.5)
  n4 = round(num_sort[2]*0.5)
  
  #within curve dependency high : time dependency high
  delta <-rep(1,time)
  delta
  t <- cumsum(delta)
  Sig2 <- Toeplitz((seq(0.3,0,by=-0.3/(time-1))))
  intercept1 <- rnorm(n1, mean=0, sd=0.1)
  intercept2 <- rnorm(n2, mean=0, sd=0.1)
  intercept3 <- rnorm(n3, mean=0, sd=0.1)
  intercept4 <- rnorm(n4, mean=0, sd=0.1)
  
  # Normal renogram curve
  normal_renogram <- function(t) {
    y <- (1 - exp(-0.6 * t)) * exp(-0.05 * (t - 10)) # 10분에서 빠르게 꺾이도록 설정
    #y <- (1 - exp(-0.5 * t)) * exp(-0.05 * (t - 5))
    return(y)
  }
  
  # Obstructive renogram curve
  obstructive_renogram <- function(t) {
    y <- (1 - exp(-0.1 * t)) + 0.01 * t # 부드럽게 증가
    return(y)
  }
  
  # Dilated unobstructive renogram curve
  dilated_unobstructive_renogram <- function(t) {
    y <- (1 - exp(-0.1 * t)) * exp(-0.03 * (t - 20)) # 20분에서 부드럽게 감소
    return(y)
  }
  
  dilated_unobstructive_2 <- function(t){
    start <- 0.7
    end <- 0.2
    range_t <- max(t) - min(t)
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.05 * sin(2 * pi * (t - 60) / range_t)
    return(y)
  }
  
  # Partial obstruction renogram curve
  partial_obstruction_renogram <- function(t) {
    y <- (1 - exp(-0.1 * t)) * exp(-0.02 * (t - 25)) # 25분에서 부드럽게 감소
    return(y)
  }
  
  #red : obstructive
  #blue : nonobstructive
  
  par(mfrow=c(2,2))
  plot(normal_renogram(c(1:time)), type="l", col="blue", main="Mean function: 1st Data", ylim=c(0,2), ylab="Value", xlab="Time")
  lines(dilated_unobstructive_renogram(c(1:time)), col="blue")
  lines(obstructive_renogram(c(1:time)), col="red")
  lines(partial_obstruction_renogram(c(1:time)), col="red")
  
  x <- rmvnorm(n=n1, mean=obstructive_renogram(t), sigma=Sig2) + intercept1 ; x2 <- rep(1,n1)
  y <- rmvnorm(n=n2, mean=partial_obstruction_renogram(t), sigma=Sig2) + intercept2  ; y2 <- rep(1,n2)
  z <- rmvnorm(n=n3, mean=dilated_unobstructive_renogram(t), sigma=Sig2) + intercept3; z2 <- rep(-1,n3)
  w <- rmvnorm(n=n4, mean=normal_renogram(t), sigma=Sig2) + intercept4 ; w2 <- rep(-1,n4)
  
  delta <-rep(1,time2)
  t <- cumsum(delta)
  Sig2 <- Toeplitz((seq(0.3,0,by=-0.3/(time2-1))))
  
  normal_renogram_2 <- function(t){
    start <- 0.2
    end <- 0.1
    range_t <- max(t) - min(t)
    mid_t <- min(t) + range_t / 2
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t))
    return(y)
  }
  
  
  obstructive_renogram2 <- function(t){
    start <- 1.5
    end <- 1.0
    range_t <- max(t) - min(t)
    mid_t <- min(t) + range_t / 2
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.05 * sin(1 * pi * (t - min(t)) / range_t)
    #y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.1 * sin(2 * pi * (t - min(t)) / range_t)
    
    return(y)
  }
  
  dilated_unobstructive_2 <- function(t){
    start <- 0.7
    end <- 0.2
    range_t <- max(t) - min(t)
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.05 * sin(2 * pi * (t - 60) / range_t)
    return(y)
  }
  
  partial_obstruction2 <-function(t){
    start <- 0.7
    end <- 0.65
    range_t <- max(t) - min(t)
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.01 * sin(2 * pi * (t - min(t)) / range_t)
    return(y)
  }
  
  
  plot(normal_renogram_2(c(1:time2)), type="l", col="blue", main="Mean function: 2nd Data", ylim=c(0,2), ylab="Value", xlab="Time")
  lines(dilated_unobstructive_2(c(1:time2)), col="blue")
  lines(obstructive_renogram2(c(1:time2)), col="red")
  lines(partial_obstruction2(c(1:time2)), col="red")
  
  x_2 <- rmvnorm(n=n1, mean=obstructive_renogram2(t), sigma=Sig2) + intercept1 ; x3 <- rep(1,n1)
  y_2 <- rmvnorm(n=n2, mean=partial_obstruction2(t), sigma=Sig2) + intercept2 ; y3 <- rep(1,n2)
  z_2 <- rmvnorm(n=n3, mean=dilated_unobstructive_2(t), sigma=Sig2) + intercept3 ; z3 <- rep(-1,n3)
  w_2 <- rmvnorm(n=n4, mean=normal_renogram_2(t), sigma=Sig2) + intercept4 ; w3 <- rep(-1,n4)
  # 
  matplot(t(z), type="l", col="blue",xlab="Time", ylab="Value", main="All Curves", ylim=c(-1,2.5))
  matplot(t(w), type="l", col="blue", add = TRUE)
  matplot(t(x), type="l", col="red", add=TRUE)
  matplot(t(y), type="l", col="red", add = TRUE)
  
  legend("topright",inset=c(0.01, 0.01),legend=c("Group 1","Group 2"),fill=c("red","blue"),cex=0.8,box.lty=0, box)
  #
  matplot(t(z_2), type="l", col="blue",xlab="Time", ylab="Value", main="All Curves", ylim=c(-1,2.5))
  matplot(t(w_2), type="l", col="blue", add = TRUE)
  matplot(t(x_2), type="l", col="red", add=TRUE)
  matplot(t(y_2), type="l", col="red", add = TRUE)
  
  legend("topright",inset=c(0.01, 0.01),legend=c("Group 1","Group 2"),fill=c("red","blue"),cex=0.8,box.lty=0, box)
  
  X <- rbind(x,y,z,w)
  X2 <- rbind(x_2,y_2,z_2,w_2)
  y <- c(x2, y2,z2,w2)
  y2 <- c(x3, y3,z3,w3)
  
  #with Y value
  XmaF <- cbind(X,y)
  XmaF2 <- cbind(X2,y2)
  
  return(list(XmaF,XmaF2))
} 
dataD<-  function(num, time, time2,iter){
  num_sort = sort(random_numbers_with_sum(num,iter))
  print(num_sort)
  n1 = round(num_sort[1]*0.5)
  n2 = round(num_sort[1]*0.5)
  n3 = round(num_sort[2]*0.5)
  n4 = round(num_sort[2]*0.5)
  
  #within curve dependency high : time dependency high
  delta <-rep(1,time)
  delta
  t <- cumsum(delta)
  Sig2 <- Toeplitz((seq(0.3,0,by=-0.3/(time-1))))
  intercept1 <- rnorm(n1, mean=0, sd=0.1)
  intercept2 <- rnorm(n2, mean=0, sd=0.1)
  intercept3 <- rnorm(n3, mean=0, sd=0.1)
  intercept4 <- rnorm(n4, mean=0, sd=0.1)
  
  # Normal renogram curve
  normal_renogram <- function(t) {
    #y <- (1 - exp(-0.5 * t)) * exp(-0.1 * (t - 5)) # 10분에서 빠르게 꺾이도록 설정
    y <- (1 - exp(-0.08 * t)) * exp(-0.01 * (t - 15))
    return(y)
  }
  
  
  # Obstructive renogram curve
  obstructive_renogram <- function(t) {
    y <- (1 - exp(-0.1 * t)) + 0.005 * t # 부드럽게 증가
    return(y)
  }
  
  
  # Dilated unobstructive renogram curve
  dilated_unobstructive_renogram <- function(t) {
    y <- (1 - exp(-0.1 * t)) * exp(-0.03 * (t - 20)) # 20분에서 부드럽게 감소
    return(y)
  }
  
  # Partial obstruction renogram curve
  partial_obstruction_renogram <- function(t) {
    y <- (1 - exp(-0.1 * t)) * exp(-0.02 * (t - 25)) # 25분에서 부드럽게 감소
    return(y)
  }
  
  par(mfrow=c(2,2))
  plot(normal_renogram(c(1:time)), type="l", col="blue", main="Mean function: 1st Data", ylim=c(0,2), ylab="Value", xlab="Time")
  lines(dilated_unobstructive_renogram(c(1:time)), col="blue")
  lines(obstructive_renogram(c(1:time)), col="red")
  lines(partial_obstruction_renogram(c(1:time)), col="red")
  
  
  x <- rmvnorm(n=n1, mean=obstructive_renogram(t), sigma=Sig2) + intercept1 ; x2 <- rep(1,n1)
  y <- rmvnorm(n=n2, mean=partial_obstruction_renogram(t), sigma=Sig2)+ intercept2 ; y2 <- rep(1,n2)
  z <- rmvnorm(n=n3, mean=dilated_unobstructive_renogram(t), sigma=Sig2) + intercept3 ; z2 <- rep(-1,n3)
  w <- rmvnorm(n=n4, mean=normal_renogram(t), sigma=Sig2) + intercept4 ; w2 <- rep(-1,n4)
  
  delta <-rep(1,time2)
  t <- cumsum(delta)
  Sig2 <- Toeplitz((seq(0.3,0,by=-0.3/(time2-1))))
  
  normal_renogram_2 <- function(t){
    start <- 0.2
    end <- 0.1
    range_t <- max(t) - min(t)
    mid_t <- min(t) + range_t / 2
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t))
    return(y)
  }
  
  obstructive_renogram2 <- function(t){
    start <- 1.3
    end <- 1.0
    range_t <- max(t) - min(t)
    mid_t <- min(t) + range_t / 2
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.1 * sin(1 * pi * (t - min(t)) / range_t)
    #y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.1 * sin(2 * pi * (t - min(t)) / range_t)
    
    return(y)
  }
  
  partial_obstruction2 <-function(t){
    start <- 0.7
    end <- 0.65
    range_t <- max(t) - min(t)
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.01 * sin(2 * pi * (t - min(t)) / range_t)
    return(y)
  }
  
  dilated_unobstructive_2 <- function(t){
    start <- 0.7
    end <- 0.2
    range_t <- max(t) - min(t)
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.05 * sin(2 * pi * (t - 60) / range_t)
    return(y)
  }
  
  
  plot(normal_renogram_2(c(1:time2)), type="l", col="blue", main="Mean function: 2nd Data", ylim=c(0,2), ylab="Value", xlab="Time")
  lines(dilated_unobstructive_2(c(1:time2)), col="blue")
  lines(obstructive_renogram2(c(1:time2)), col="red")
  lines(partial_obstruction2(c(1:time2)), col="red")
  
  x_2 <- rmvnorm(n=n1, mean=obstructive_renogram2(t), sigma=Sig2) + intercept1; x3 <- rep(1,n1)
  y_2 <- rmvnorm(n=n2, mean=partial_obstruction2(t), sigma=Sig2) + intercept2 ; y3 <- rep(1,n2)
  z_2 <- rmvnorm(n=n3, mean=dilated_unobstructive_2(t), sigma=Sig2) + intercept3 ; z3 <- rep(-1,n3)
  w_2 <- rmvnorm(n=n4, mean=normal_renogram_2(t), sigma=Sig2) + intercept4 ; w3 <- rep(-1,n4)
  
  matplot(t(z), type="l", col="blue",xlab="Time", ylab="Value", main="All Curves", ylim=c(0,2))
  matplot(t(w), type="l", col="blue", add = TRUE)
  matplot(t(x), type="l", col="red", add=TRUE)
  matplot(t(y), type="l", col="red", add = TRUE)
  
  legend("bottomright",legend=c("Obstruct","Normal"),fill=c("red","blue"),cex=0.8)
  
  matplot(t(z_2), type="l", col="blue",xlab="Time", ylab="Value", main="All Curves", ylim=c(0,2))
  matplot(t(w_2), type="l", col="blue", add = TRUE)
  matplot(t(x_2), type="l", col="red", add=TRUE)
  matplot(t(y_2), type="l", col="red", add = TRUE)
  
  legend("bottomright",legend=c("Obstruct","Normal"),fill=c("red","blue"),cex=0.8)
  
  X <- rbind(x,y,z,w)
  X2 <- rbind(x_2,y_2,z_2,w_2)
  y <- c(x2, y2,z2,w2)
  y2 <- c(x3, y3,z3,w3)
  
  #with Y value
  XmaF <- cbind(X,y)
  XmaF2 <- cbind(X2,y2)
  
  return(list(XmaF,XmaF2))
}
dataE <- function(num, time, time2,iter){
  set.seed(iter)
  num_sort = sort(random_numbers_with_sum(num,iter))
  print(num_sort)
  n1 = round(num_sort[1]*0.5)
  n2 = round(num_sort[1]*0.5)
  n3 = round(num_sort[2]*0.5)
  n4 = round(num_sort[2]*0.5)
  
  #within curve dependency high : time dependency high
  delta <-rep(1,time)
  delta
  t <- cumsum(delta)
  Sig2 <- Toeplitz((seq(0.3,0,by=-0.3/(time-1))))
  intercept1 <- rnorm(n1, mean=0, sd=0.1)
  intercept2 <- rnorm(n2, mean=0, sd=0.1)
  intercept3 <- rnorm(n3, mean=0, sd=0.1)
  intercept4 <- rnorm(n4, mean=0, sd=0.1)
  
  # Normal renogram curve
  normal_renogram <- function(t) {
    return(3*(0.3*sin((pi/18)* t- pi/4)+ 0.5))
  }
  
  # Obstructive renogram curve
  obstructive_renogram <- function(t) {
    return(3*(0.2 * sin((pi/17) * t - pi/4) + 0.5))
  }
  
  # Dilated unobstructive renogram curve
  dilated_unobstructive_renogram <- function(t) {
    return(3*(0.2*sin((pi/15)* t- pi/4)+ 0.5))
    
  }
  
  # Partial obstruction renogram curve
  partial_obstruction_renogram <- function(t) {
    return(3*(0.3*sin((pi/20)* t- pi/4)+ 0.5))
    
  }
  
  
  par(mfrow=c(2,2))
  plot(normal_renogram(c(1:time)), type="l", col="blue", main="Mean function: 1st Data", ylab="Value", xlab="Time", ylim=c(0,4))
  lines(dilated_unobstructive_renogram(c(1:time)), col="blue")
  lines(obstructive_renogram(c(1:time)), col="red")
  lines(partial_obstruction_renogram(c(1:time)), col="red")
  
  x <- rmvnorm(n=n1, mean=obstructive_renogram(t), sigma=Sig2) + intercept1 ; x2 <- rep(-1,n3)
  y <- rmvnorm(n=n2, mean=partial_obstruction_renogram(t), sigma=Sig2) + intercept2  ; y2 <- rep(-1,n4)
  z <- rmvnorm(n=n3, mean=dilated_unobstructive_renogram(t), sigma=Sig2) + intercept3; z2 <- rep(1,n1)
  w <- rmvnorm(n=n4, mean=normal_renogram(t), sigma=Sig2) + intercept4 ; w2 <- rep(1,n2)
  
  delta <-rep(1,time2)
  t <- cumsum(delta)
  Sig2 <- Toeplitz((seq(0.3,0,by=-0.3/(time2-1))))
  
  normal_renogram_2 <- function(t){
    start <- 0.8
    end <- 0.4
    range_t <- max(t) - min(t)
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.01 * sin(2 * pi * (t - min(t)) / range_t)
    return(3*y)      
  }
  
  
  obstructive_renogram2 <- function(t){
    start <- 0.3
    end <- 0.2
    range_t <- max(t) - min(t)
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.01 * sin(2 * pi * (t - min(t)) / range_t)
    return(3*y)    
  }
  
  dilated_unobstructive_2 <- function(t){
    #0.2 * sin((pi / 70) * t) + 0.6
    start <- 0.5
    end <- 0.4
    range_t <- max(t) - min(t)
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.01 * sin(2 * pi * (t - min(t)) / range_t)
    return(3*y)
  }
  
  partial_obstruction2 <-function(t){
    start <- 0.6
    end <- 0.3
    range_t <- max(t) - min(t)
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.01 * sin(2 * pi * (t - min(t)) / range_t)
    return(3*y)
  }
  
  
  plot(normal_renogram_2(c(1:time2)), type="l", col="blue" ,main="Mean function: 2nd Data", ylab="Value", xlab="Time", ylim=c(0,4))
  lines(dilated_unobstructive_2(c(1:time2)), col="blue")
  lines(obstructive_renogram2(c(1:time2)), col="red")
  lines(partial_obstruction2(c(1:time2)), col="red")
  
  x_2 <- rmvnorm(n=n1, mean=obstructive_renogram2(t), sigma=Sig2) + intercept1 ; x3 <- rep(-1,n3)
  y_2 <- rmvnorm(n=n2, mean=partial_obstruction2(t), sigma=Sig2) + intercept2 ; y3 <- rep(-1,n4)
  z_2 <- rmvnorm(n=n3, mean=dilated_unobstructive_2(t), sigma=Sig2) + intercept3 ; z3 <- rep(1,n1)
  w_2 <- rmvnorm(n=n4, mean=normal_renogram_2(t), sigma=Sig2) + intercept4 ; w3 <- rep(1,n2)
  # 
  matplot(t(z), type="l", col="blue",xlab="Time", ylab="Value", main="All Curves")
  matplot(t(w), type="l", col="blue", add = TRUE)
  matplot(t(x), type="l", col="red", add=TRUE)
  matplot(t(y), type="l", col="red", add = TRUE)
  
  legend("topright",inset=c(0.01, 0.01),legend=c("Group 1","Group 2"),fill=c("red","blue"),cex=0.8,box.lty=0, box)
  #
  matplot(t(z_2), type="l", col="blue",xlab="Time", ylab="Value", main="All Curves")
  matplot(t(w_2), type="l", col="blue", add = TRUE)
  matplot(t(x_2), type="l", col="red", add=TRUE)
  matplot(t(y_2), type="l", col="red", add = TRUE)
  
  legend("topright",inset=c(0.01, 0.01),legend=c("Group 1","Group 2"),fill=c("red","blue"),cex=0.8,box.lty=0, box)
  
  X <- rbind(x,y,z,w)
  X2 <- rbind(x_2,y_2,z_2,w_2)
  y <- c(x2, y2,z2,w2)
  y2 <- c(x3, y3,z3,w3)
  XmaF <- cbind(X,y)
  XmaF2 <- cbind(X2,y2)
  
  return(list(XmaF,XmaF2))
}

dataF <- function(num, time, time2,iter){
  set.seed(iter)
  num_sort = sort(random_numbers_with_sum(num,iter))
  print(num_sort)
  n1 = round(num_sort[1]*0.5)
  n2 = round(num_sort[1]*0.5)
  n3 = round(num_sort[2]*0.5)
  n4 = round(num_sort[2]*0.5)
  intercept1 <- rnorm(n1, mean=0, sd=0.1)
  intercept2 <- rnorm(n2, mean=0, sd=0.1)
  intercept3 <- rnorm(n3, mean=0, sd=0.1)
  intercept4 <- rnorm(n4, mean=0, sd=0.1)
  delta <-rep(1,time)
  delta
  t <- cumsum(delta)
  Sig2 <- Toeplitz((seq(0.3,0,by=-0.3/(time-1))))
  
  # Normal renogram curve
  normal_renogram <- function(t) {
    return(3*(0.3*sin((pi/18)* t- pi/4)+ 0.5))
  }
  
  # Obstructive renogram curve
  obstructive_renogram <- function(t) {
    return(3*(0.2 * sin((pi/17) * t - pi/4) + 0.5))
  }
  
  # Dilated unobstructive renogram curve
  dilated_unobstructive_renogram <- function(t) {
    return(3*(0.2*sin((pi/15)* t- pi/4)+ 0.5))
    
  }
  
  # Partial obstruction renogram curve
  partial_obstruction_renogram <- function(t) {
    return(3*(0.3*sin((pi/20)* t- pi/4)+ 0.5))
    
  }
  
  
  par(mfrow=c(2,2))
  plot(normal_renogram(c(1:time)), type="l", col="blue", main="Mean function: 1st Data", ylab="Value", xlab="Time", ylim=c(0,4))
  lines(dilated_unobstructive_renogram(c(1:time)), col="blue")
  lines(obstructive_renogram(c(1:time)), col="red")
  lines(partial_obstruction_renogram(c(1:time)), col="red")
  
  x <- rmvnorm(n=n1, mean=obstructive_renogram(t), sigma=Sig2) +intercept1; x2 <- rep(-1,n3)
  y <- rmvnorm(n=n2, mean=partial_obstruction_renogram(t), sigma=Sig2) + intercept2  ; y2 <- rep(-1,n4)
  z <- rmvnorm(n=n3, mean=dilated_unobstructive_renogram(t), sigma=Sig2) + intercept3 ; z2 <- rep(1,n1)
  w <- rmvnorm(n=n4, mean=normal_renogram(t), sigma=Sig2)+ intercept4 ; w2 <- rep(1,n2)
  
  delta <-rep(1,time2)
  t <- cumsum(delta)
  Sig2 <- Toeplitz((seq(0.3,0,by=-0.3/(time2-1))))
  
  normal_renogram_2 <- function(t){
    #0.2 * sin((pi / 30) * x) + 0.75
    start <- 0.8
    end <- 0.4
    range_t <- max(t) - min(t)
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.02 * sin(2 * pi * (t - min(t)) / range_t)
    return(3*y)    
  }
  
  
  obstructive_renogram2 <- function(t){
    start <- 0.4
    end <- 0.2
    range_t <- max(t) - min(t)
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.01 * sin(2 * pi * (t - min(t)) / range_t)
    return(3*y)    
  }
  
  dilated_unobstructive_2 <- function(t){
    start <- 0.6
    end <- 0.4
    range_t <- max(t) - min(t)
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.01 * sin(2 * pi * (t - min(t)) / range_t)
    return(3*y)
  }
  
  partial_obstruction2 <-function(t){
    3*(0.3 * sin((pi / 40) * t) + 0.5)
  }
  
  
  plot(normal_renogram_2(c(1:time2)), type="l", col="blue",main="Mean function: 2nd Data", ylab="Value", xlab="Time", ylim=c(0,4))
  lines(dilated_unobstructive_2(c(1:time2)), col="blue")
  lines(obstructive_renogram2(c(1:time2)), col="red")
  lines(partial_obstruction2(c(1:time2)), col="red")
  
  x_2 <- rmvnorm(n=n1, mean=obstructive_renogram2(t), sigma=Sig2) + intercept1 ; x3 <- rep(-1,n3)
  y_2 <- rmvnorm(n=n2, mean=partial_obstruction2(t), sigma=Sig2) + intercept2 ; y3 <- rep(-1,n4)
  z_2 <- rmvnorm(n=n3, mean=dilated_unobstructive_2(t), sigma=Sig2) + intercept3; z3 <- rep(1,n1)
  w_2 <- rmvnorm(n=n4, mean=normal_renogram_2(t), sigma=Sig2) + intercept4 ; w3 <- rep(1,n2)
  # 
  matplot(t(z), type="l", col="blue",xlab="Time", ylab="Value", main="All Curves")
  matplot(t(w), type="l", col="blue", add = TRUE)
  matplot(t(x), type="l", col="red", add=TRUE)
  matplot(t(y), type="l", col="red", add = TRUE)
  
  legend("topright",inset=c(0.01, 0.01),legend=c("Group 1","Group 2"),fill=c("red","blue"),cex=0.8,box.lty=0, box)
  #
  matplot(t(z_2), type="l", col="blue",xlab="Time", ylab="Value", main="All Curves")
  matplot(t(w_2), type="l", col="blue", add = TRUE)
  matplot(t(x_2), type="l", col="red", add=TRUE)
  matplot(t(y_2), type="l", col="red", add = TRUE)
  
  legend("topright",inset=c(0.01, 0.01),legend=c("Group 1","Group 2"),fill=c("red","blue"),cex=0.8,box.lty=0, box)
  
  X <- rbind(x,y,z,w)
  X2 <- rbind(x_2,y_2,z_2,w_2)
  y <- c(x2, y2,z2,w2)
  y2 <- c(x3, y3,z3,w3)
  
  XmaF <- cbind(X,y)
  XmaF2 <- cbind(X2,y2)
  
  return(list(XmaF,XmaF2))
}
dataB <- function(num, time, time2,iter){
  set.seed(iter)

  num_sort = sort(random_numbers_with_sum(num,iter))
  print(num_sort)
  n1 = round(num_sort[1])
  n2 = round(num_sort[2])
  
  delta <-rep(1,time)
  delta
  t <- cumsum(delta)
  Sig2 <- Toeplitz((seq(0.3,0,by=-0.3/(time-1))))
  intercept1 <- rnorm(n1, mean=0, sd=0.1)
  intercept2 <- rnorm(n2, mean=0, sd=0.1)
  
  curve1 <- function(x) {
    return(3*(0.1 * sin((pi/35) * x - pi/2) + 0.5))
  }
  
  curve2 <- function(x){
    return(3*(0.1*sin((pi/40)* x- pi/4)+ 0.5))
  }

  x <- rmvnorm(n=n2, mean=curve1(t), sigma=Sig2) + intercept2 ; x2 <- rep(-1,n2)
  y <- rmvnorm(n=n1, mean=curve2(t), sigma=Sig2) + intercept1 ; y2 <- rep(1,n1)
  
  delta <-rep(1,time2)
  t <- cumsum(delta)
  Sig2 <- Toeplitz((seq(0.3,0,by=-0.3/(time2-1))))
  
  curve3 <- function(x) {
    start <- 0.7
    end <- 0.3
    range_t <- max(t) - min(t)
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.02 * sin(2 * pi * (t - min(t)) / range_t)
    return(3*y)
  }
  
  curve4 <- function(x) {
    return((3*(0.05*sin((pi/40)*(x-40)-pi/3)+ 0.5)))
    
  }
  
  x_2 <- rmvnorm(n=n2, mean=curve3(t), sigma=Sig2) + intercept2 ; x3 <- rep(-1,n2)
  y_2 <- rmvnorm(n=n1, mean=curve4(t), sigma=Sig2) + intercept1 ; y3 <- rep(1,n1)
  
  par(mfrow=c(2,2))
  plot(curve1(c(1:time)), type="l", col="blue", main="Mean function: 1st Data",ylab="Value", xlab="Time")
  lines(curve2(c(1:time)), col="red",ylim=c(0,1))
  
  plot(curve3(c(1:time2)), type="l", col="blue", main="Mean function: 2nd Data",ylab="Value", xlab="Time")
  lines(curve4(c(1:time2)), col="red", ylim=c(0,1))
  
  
  matplot(t(x), type="l", col="blue",xlab="Time", ylab="Value", main="All Curves")
  matplot(t(y), type="l", col="red", add = TRUE, ylim=c(-1,2))
  
  legend("topright", inset=c(0.01, 0.01), legend=c("Group 1","Group 2"),fill=c("red","blue"),cex=0.8, box.lty=0)
  
  matplot(t(x_2), type="l", col="blue",xlab="Time", ylab="Value", main="All Curves")
  matplot(t(y_2), type="l", col="red", add = TRUE, ylim=c(-1,2))
  
  
  legend("topright",inset=c(0.01, 0.01),legend=c("Group 1","Group 2"),fill=c("red","blue"),cex=0.8,box.lty=0, box)

  X <- rbind(x,y)
  X2 <- rbind(x_2,y_2)
  y <- c(x2, y2)
  y2 <- c(x3, y3)
  
  XmaF <- cbind(X,y)
  XmaF2 <- cbind(X2,y2)
  
  return(list(XmaF,XmaF2))
}

dataA <- function(num, time, time2,iter){

  set.seed(iter)
  
  num_sort = sort(random_numbers_with_sum(num,iter))
  print(num_sort)
  n1 = round(num_sort[1])
  n2 = round(num_sort[2])
  
  #within curve dependency high : time dependency high
  delta <-rep(1,time)
  delta
  t <- cumsum(delta)
  Sig2 <- Toeplitz((seq(0.3,0,by=-0.3/(time-1))))
  intercept1 <- rnorm(n1, mean=0, sd=0.1)
  intercept2 <- rnorm(n2, mean=0, sd=0.1)
  
  curve1 <- function(x) {
    return((sin((pi/35) * x - pi/2) + 0.5))
  }
  
  curve2 <- function(x){
    return((sin((pi/35) * x - pi/4) + 0.5))
  }
  
  x <- rmvnorm(n=n1, mean=curve1(t), sigma=Sig2) + intercept1 ; x2 <- rep(1,n1)
  y <- rmvnorm(n=n2, mean=curve2(t), sigma=Sig2) + intercept2 ; y2 <- rep(-1,n2)
  
  delta <-rep(1,time2)
  t <- cumsum(delta)
  Sig2 <- Toeplitz((seq(0.3,0,by=-0.3/(time2-1))))
  
  curve3 <- function(x) {
    start <- 1.5
    end <- 0
    range_t <- max(t) - min(t)
    amplitude <- (start - end) / 2
    y <- end + amplitude * (1 + cos(pi * (t - min(t)) / range_t)) + 0.02 * sin(2 * pi * (t - min(t)) / range_t)
    return(y)
  }
  
  curve4 <- function(x) {
    return(((0.1*sin((pi/40)*(x-30)- pi/4)+ 0.5)))
  }
  
  x_2 <- rmvnorm(n=n1, mean=curve3(t), sigma=Sig2) + intercept1 ; x3 <- rep(1,n1)
  y_2 <- rmvnorm(n=n2, mean=curve4(t), sigma=Sig2) + intercept2 ; y3 <- rep(-1,n2)
  
  par(mfrow=c(2,2))
  plot(curve1(c(1:time)), type="l", col="blue", main="Mean function: 1st Data",ylab="Value", xlab="Time", ylim=c(-1,3))
  lines(curve2(c(1:time)), col="red")
  
  plot(curve3(c(1:time2)), type="l", col="blue", main="Mean function: 2nd Data",ylab="Value", xlab="Time", ylim=c(-1,3))
  lines(curve4(c(1:time2)), col="red", ylim=c(0,1))
  
  
  matplot(t(x), type="l", col="blue",xlab="Time", ylab="Value", main="All Curves")
  matplot(t(y), type="l", col="red", add = TRUE, ylim=c(-1,2))
  
  legend("topright", inset=c(0.01, 0.01), legend=c("Group 1","Group 2"),fill=c("red","blue"),cex=0.8, box.lty=0)
  
  matplot(t(x_2), type="l", col="blue",xlab="Time", ylab="Value", main="All Curves")
  matplot(t(y_2), type="l", col="red", add = TRUE, ylim=c(-1,2))
  
  
  legend("topright",inset=c(0.01, 0.01),legend=c("Group 1","Group 2"),fill=c("red","blue"),cex=0.8,box.lty=0, box)
  
  X <- rbind(x,y)
  X2 <- rbind(x_2,y_2)
  y <- c(x2, y2)
  y2 <- c(x3, y3)
  
 XmaF <- cbind(X,y)
  XmaF2 <- cbind(X2,y2)
  
  return(list(XmaF,XmaF2))
}
