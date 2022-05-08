#programmation dynamique
opti_seg <- function(y, C, K){
  L <- length(y)
  Q0 <- matrix(0, ncol = L, nrow = L)
  for (j in 1:L){
    for (i in 1:j){
      Q0[i,j]<-C(y,i,j)
    }
  }
  Q <- matrix(0,ncol = L, nrow = K+1)
  Q[1,] <- Q0[1,]
  tau <- matrix(0,ncol = K, nrow = K)
  v_tau <- (1:(L-1))
  for (k in 1:K){
    for (j in 2:L){
      v <- 1:(j-1)
      Q[(k+1),j] <- min(Q[k,v]+Q0[(v+1),j])
      tau[k,1] <- which.min( Q[k,v_tau] + Q0[(v_tau+1),L] )
    }
    if (k>=2){
      for (i in 2:k){
        tau[k,i] <- which.min(Q[(k-i+1),v_tau] + Q0[(v_tau+1),tau[k,(i-1)]])
      }
    }
  }
  return(list("tau"=tau,"couts_finaux"=Q[,L][-1],"Q"=Q,"Q0"=Q0))
}

opti_seg_ar1 <- function(y, C, K){
  L <- length(y)
  Q0 <- matrix(0, ncol = L-1, nrow = L-1)
  for (j in 1:(L-1)){
    for (i in 1:j){
      Q0[i,j]<-C(y,i+1,j+1)
    }
  }
  Q <- matrix(0,ncol = L-1, nrow = K+1)
  Q[1,] <- Q0[1,]
  tau <- matrix(0,ncol = K, nrow = K)
  v_tau <- (2:(L-1))
  for (k in 1:K){
    for (j in 2:(L-1)){
      v <- 2:j
      Q[k+1,j] <- min(Q[k,v-1]+Q0[v,j])
      tau[k,1] <- which.min( Q[k,v_tau-1] + Q0[v_tau,L-1] )
    }
    if (k>=2){
      for (i in 2:k){
        tau[k,i] <- which.min( Q[(k-i+1),v_tau-1] + Q0[(v_tau),tau[k,(i-1)]-1] )
      }
    }
  }
  return(list("tau"=tau,"couts_finaux"=Q[-1,(L-1)],"Q"=Q,"Q0"=Q0))
}

#PELT
pelt_mean <- function(y, min.seg.len, pen, cout)
{
  n <- length(y)
  FF <- rep(0, n)
  pt_rupt  <- rep(0, n)  
  pt_rupt[1:(2*min.seg.len)] <- 0
  FF[1] <- -pen
  FF[2:(2*min.seg.len)] <-sapply(2:(2*min.seg.len),function(x){cout(y,1,x)})
  R <- c(1,min.seg.len)
  for (s in (2*min.seg.len):n){
    R.inter <- R[R<=(s-min.seg.len)]
    cout.inter <-sapply(R.inter,function(x){cout(y,x,s)})
    templike <- FF[R.inter]+ cout.inter + pen
    FF[s] <- min(templike)
    pt_rupt[s] <- R[which.min(templike)]
    R <- c(R[templike-pen < FF[s]], s)
  }
  i <- n
  pt_rup_final <- c()
  while(i>1){
    if (pt_rupt[i]>1){
      pt_rup_final <- c(pt_rup_final, pt_rupt[i])
    }
    i <- pt_rupt[i]-1
  }
  return(pt_rup_final)
}

pelt_trend <- function(y, min.seg.len, pen, cout)
{
  n <- length(y)
  FF <- rep(0, n)
  pt_rupt  <- rep(0, n)  
  pt_rupt[1:(2*min.seg.len)] <- 0
  FF[1] <- -pen
  FF[2:(2*min.seg.len)] <-sapply(2:(2*min.seg.len),function(x){cout(y,1,x)})
  R <- c(1,min.seg.len)
  for (s in (2*min.seg.len):n){
    R.inter <- R[R<=(s-min.seg.len)]
    cout.inter <-sapply(R.inter,function(x){cout(y,x,s)})
    templike <- FF[R.inter]+ cout.inter + pen
    FF[s] <- min(templike)
    pt_rupt[s] <- R[which.min(templike)]
    R <- c(R[templike-pen < FF[s]], s)
  }
  i <- n
  pt_rup_final <- c()
  while(i>1){
    if (pt_rupt[i]>1){
      pt_rup_final <- c(pt_rup_final, pt_rupt[i])
    }
    i <- pt_rupt[i]-1
  }
  return(pt_rup_final)
}

pelt_AR <- function(y, min.seg.len, pen, cout)
  
{
  n <- length(y)
  FF <- rep(0, n)
  pt_rupt  <- rep(0, n)  
  pt_rupt[1:(2*min.seg.len)] <- 0
  FF[1] <- -pen
  FF[2:(2*min.seg.len)] <-sapply(2:(2*min.seg.len),function(x){cout(y,2,x)})
  R <- c(2,min.seg.len)
  for (s in (2*min.seg.len):n){
    R.inter <- R[R<=(s-min.seg.len)]
    cout.inter <-sapply(R.inter,function(x){cout(y,x,s)})
    templike <- FF[R.inter]+ cout.inter + pen
    FF[s] <- min(templike)
    pt_rupt[s] <- R[which.min(templike)]
    R <- c(R[templike-pen < FF[s]], s)
  }
  i <- n
  pt_rup_fin <- c()
  while(i>1){
    if (pt_rupt[i]>2){
      pt_rup_fin <- c(pt_rup_fin, pt_rupt[i])
    }
    i <- pt_rupt[i]-1
  }
  return(pt_rup_fin)
}


#################################################################################
#modele 1
est_mod1 <- function(y,t1,t2){
  return(list("mu"=mean(y[t1:t2])))
}

cout_mod1 <- function(y,t1,t2){
  return(sum((y[t1:t2]-mean(y[t1:t2]))^2))
}

simu_mod1 <- function(l,m,tau){
  tau <- c(0,tau,l)
  y <- numeric(l)
  for (k in 1:(length(tau)-1)){
    epsilon <- rnorm(tau[k+1]-tau[k])
    y[(tau[k]+1):tau[k+1]] <- m[k] + epsilon
  }
  return(y)  
}

#modele 2
est_mod2 <- function(y,t1,t2){
  yt <- y[t1:bt2]
  beta <- sum((yt - mean(yt)) * ((t1:t2) - (t1 + t2 + 1)/2)) / sum(((t1:t2) - (t1 + t2 + 1)/2)^2)
  mu <- mean(yt) - beta * (t1 + t2 + 1)/2
  return(list("mu"=mu,"beta"=beta))
}

cout_mod2 <- function(y, t1, t2){
  yt <- y[t1:t2]
  beta <- sum((yt - mean(yt)) * ((t1:t2) - (t1 + t2 + 1)/2)) / sum(((t1:t2) - (t1 + t2 + 1)/2)^2)
  mu <- mean(yt) - beta * (t1 + t2 + 1)/2
  return(sum((yt - (mu + beta * t1:t2))^2))
}

simu_mod2 <- function(l,m,b,tau){
  tau <- c(0,tau,l)
  y <- numeric(l)
  for (k in 1:(length(tau)-1)){
    epsilon <- rnorm(tau[k+1]-tau[k])
    y[(tau[k]+1):tau[k+1]] <- m[k] * ((tau[k]+1):tau[k+1]) + b[k] + epsilon
  }
  return(y)
}

#simulation AR1
simu_ar1 <- function(l,p,v_){
  y0 <- 0
  for (i in 1:50) {
    y0 <- y0*p + rnorm(1,0,sqrt(v_))
  }
  y <- c(y0, numeric(l-1))
  for (t in 2:l){
    y[t] <- y[t-1]*p + rnorm(1,0,sqrt(v_))
  }
  return(y)
}

#modele 3

est_mod3 <- function(y,t1,t2){
  yt <- y[t1:t2]
  ytm1 <- y[(t1-1):(t2-1)]
  tt <- t1:t2
  I <- t2 - t1 + 1
  phi <- (I * sum(ytm1 * yt) - sum(ytm1) * sum(yt))/(I * sum(ytm1^2) - sum(ytm1)^2)
  mu <- sum(yt - phi * ytm1)/(I*(1 - phi))
  return(list("phi"=phi,"mu"=mu))
}

cout_mod3 <- function(y,t1,t2){
  if (t1 == t2){
    return(0)
  }
  yt <- y[t1:t2]
  ytm1 <- y[(t1-1):(t2-1)]
  tt <- t1:t2
  I <- t2 - t1 + 1
  b <- I * sum(ytm1 * yt) - sum(ytm1) * sum(yt)
  e <- I * sum(ytm1^2) - sum(ytm1)^2
  phi <- b/e
  mu <- sum(yt - phi * ytm1)/(I*(1 - phi))
  cout <- sum( (yt - mu - phi * (ytm1 - mu))^2 )
  return(cout)
}

simu_mod3 <- function(l,m,p,v,tau){
  simu_ar1 <- function(l,p,v_){
    y0 <- 0
    for (i in 1:50) {
      y0 <- y0*p + rnorm(1,0,sqrt(v_))
    }
    y <- c(y0, numeric(l-1))
    for (t in 2:l){
      y[t] <- y[t-1]*p + rnorm(1,0,sqrt(v_))
    }
    return(y)
  }
  tau <- c(0,tau,l)
  y <- numeric(l)
  for (k in 1:(length(tau)-1)){
    y[(tau[k]+1):tau[k+1]] <- simu_ar1 (tau[k+1]-tau[k], p[k], v[k]) + m[k]
  }
  return(y)
}

#modele 4
est_mod4 <- function(y, t1, t2){
  yt <- y[t1:t2]
  ytm1 <- y[(t1-1):(t2-1)]
  tt <- t1:t2
  I <- t2 - t1 + 1
  
  a <- I * sum(tt * ytm1) - sum(tt) * sum(ytm1)
  c <- sum(tt)^2 - I * sum(tt^2)
  d <- I * sum(tt * yt) - sum(tt) * sum(yt) - a
  
  f <- sum(tt^2) * sum(ytm1^2) - sum(tt * ytm1)^2
  g <- sum(tt^2) * sum(ytm1 * yt) - sum(tt * yt) * sum(ytm1 * tt)
  h <- sum(tt^2) * sum(ytm1) - sum(tt) * sum(ytm1 * tt)
  i <- sum(tt)^2 - I * sum(tt^2)
  j <- sum(tt^2) * sum(yt) - sum(tt) * sum(yt * tt)
  k <- sum(tt) * sum(tt * ytm1) - sum(tt^2) * sum(ytm1)
  
  
  phi <- (h*j + g*i)/(f*i - h * k)
  w <- ((d + a) - phi * a)/c
  alpha <- w/(phi - 1)
  delta <- -sum(yt - phi * ytm1 + w * tt)/I
  beta <- -(delta + alpha * phi)/(1-phi)
  
  return(list("alpha"=alpha, "beta" = beta, "phi" = phi))
  
}

cout_mod4 <- function(y, t1, t2) {
  if (t1 == t2){
    return(0)
  }
  yt <- y[t1:t2]
  ytm1 <- y[(t1-1):(t2-1)]
  tt <- t1:t2
  I <- t2 - t1 + 1
  i <- sum(tt)^2 - I * sum(tt^2)
  h <- sum(tt^2) * sum(ytm1) - sum(tt) * sum(ytm1 * tt)
  g <- sum(tt^2) * sum(ytm1 * yt) - sum(tt * yt) * sum(ytm1 * tt)
  j <- sum(tt^2) * sum(yt) - sum(tt) * sum(yt * tt)
  k <- sum(tt) * sum(tt * ytm1) - sum(tt^2) * sum(ytm1)
  f <- sum(tt^2) * sum(ytm1^2) - sum(tt * ytm1)^2
  a <- I * sum(tt * ytm1) - sum(tt) * sum(ytm1)
  c <- sum(tt)^2 - I * sum(tt^2)
  d <- I * sum(tt * yt) - sum(tt) * sum(yt) - a
  
  
  phi <- (h*j + g*i)/(f*i - h * k)
  w <- ((d + a) - phi * a)/c
  alpha <- w/(phi - 1)
  delta <- -sum(yt - phi * ytm1 + w * tt)/I
  beta <- (delta + alpha * phi)/(phi - 1)
  
  cout <- sum( (yt - phi * ytm1 - alpha * (1 - phi) * tt - alpha * phi - beta * (1 - phi))^2  )
  if(is.nan(cout)){
    return(0)
  }
  return(cout)
}

simu_mod4 <- function(l,a,b,p,v,tau){
  tau <- c(0,tau,l)
  y <- numeric(l)
  for (k in 1:(length(tau)-1)){
    y[(tau[k]+1):tau[k+1]] <- a[k] * ((tau[k]+1):tau[k+1]) + b[k] + simu_ar1(tau[k+1]-tau[k],p[k],v[k])
  }
  return(y)
}
