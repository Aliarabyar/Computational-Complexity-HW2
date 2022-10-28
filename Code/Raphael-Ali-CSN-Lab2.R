#=======================================================================
# CSN - Lab Assignment 2
# October 2022
#=======================================================================
# Raphaël Vignon
# Ali Arabyarmohammadi
#=======================================================================
# Note: Please set the working directory to the current folder

  
require("stats4") # for MLE
require("VGAM") # for the Riemann-zeta function

write_summary <- function(language,degree_sequence) {
  cat(language,length(degree_sequence$V1),max(degree_sequence$V1),sum(degree_sequence$V1)/length(degree_sequence$V1),length(degree_sequence$V1)/sum(degree_sequence$V1),"\n")
}

get_AIC <- function(m2logL,K,N) {
  m2logL + 2*K*N/(N-K-1) # AIC with a correction for sample size
}

H <- function(gamma,k_max){
  sum(seq(k_max)^-gamma)
}

zeta_dist_truncated <- function(k,gamma,k_max){
  k^(-gamma)/H(gamma,k_max)
}

zeta_dist <- function(k,gamma){
  k^(-gamma)/zeta(gamma)
}

poisson_dist <- function(k,lambda){
  lambda^k*exp(-lambda)/(factorial(k)*(1-exp(-lambda)))
} 


geometric_dist <- function(k, p) {
  dgeom(k,p)
} 

altm_dist <- function(k,p,q){
  N=length(k)
  return(((sum((1:N)^(-p)*exp(-q*(1:N))))^(-1))*k^(-p)*exp(-q*k))
}

######################### Zeta likelihood ##################



loglikelyhood_zeta <- function(degree_sequence){
  x <- degree_sequence$V1
  
  minus_log_likelihood_zeta <- function(gamma) {
    length(x) * log(zeta(gamma)) + gamma * sum(log(x))
  }
  
  mle_zeta <- mle(minus_log_likelihood_zeta,
                  start = list(gamma = 2),
                  method = "L-BFGS-B",
                  lower = c(1.0000001))
  
  list(
    gamma = attributes(summary(mle_zeta))$coef[1],
    AIC = get_AIC(attributes(summary(mle_zeta))$m2logL, 1, length(degree_sequence$V1)),
    distribution = function(x) { zeta_dist(x, attributes(summary(mle_zeta))$coef[1])}
  )
}

loglikelyhood_zeta_truncated <- function(degree_sequence){
  
  x <- degree_sequence$V1
  minus_log_likelihood_zeta_truncated <- function(gamma, k_max) {
    length(x) * log(H(gamma,k_max)) + gamma * sum(log(x))
  }
  
  
  N <- length(x)
  MaxDeg=max(degree_sequence$V1)
  
  mle_zeta_truncated <- mle(minus_log_likelihood_zeta_truncated,
                            start = list(gamma = 2, k_max = MaxDeg),
                            method = "L-BFGS-B",
                            lower = c(1.0000001, N))
  
  summary(mle_zeta_truncated)
  
  list(
    gamma = attributes(summary(mle_zeta_truncated))$coef[1],
    kmax = attributes(summary(mle_zeta_truncated))$coef[2],
    AIC = get_AIC(attributes(summary(mle_zeta_truncated))$m2logL, 2, length(degree_sequence$V1)),
    distribution = function(x) { zeta_dist_truncated(x, attributes(summary(mle_zeta_truncated))$coef[1], attributes(summary(mle_zeta_truncated))$coef[2])}
  )
}

loglikelyhood_zeta_gamma2 <- function(degree_sequence){
  x <- degree_sequence$V1
  
  minus_log_likelihood_zeta <- function(gamma) {
    length(x) * log(zeta(gamma)) + gamma * sum(log(x))
  }
  m2logL <- 2 * minus_log_likelihood_zeta(2)
  list(
    AIC = get_AIC(m2logL, 0, length(degree_sequence$V1)),
    distribution = function(x) { zeta_dist(x, 2)}
  )
}

######################### Poisson likelihood ##################

loglikelyhood_poisson <- function(degree_sequence){
  x <- degree_sequence$V1
  N <- length(x)
  M <- sum(x)
  C <- sum(sapply(x,function(j) sum(log(2:j))))
  
  minus_log_likelihood_poisson <- function(lambda) {
    N * (lambda + log(1 - exp(-lambda))) - M * log(lambda) + C
  }
  
  mle_poisson <- mle( minus_log_likelihood_poisson,
                      start = list(lambda = M / N),
                      method = "L-BFGS-B",
                      lower = c(1.0000001) )
  
  list(
    k = attributes(summary(mle_poisson))$coef[1],
    AIC = get_AIC(attributes(summary(mle_poisson))$m2logL, 1, N),
    distribution = function(x) {poisson_dist(x, attributes(summary(mle_poisson))$coef[1])}
  )
  
}


######################### Geometric loglikelyhood ##################

loglikelyhood_geometric <- function(degree_sequence){
  x <- degree_sequence$V1
  N <- length(x)
  M <- sum(x)
  C <- sum(sapply(x,function(j) sum(log(2:j))))
  
  minus_log_likelihood_geometric <- function(q) {
    (N - M) * log(1 - q) - N * log(q)
  }
  
  mle_geometric <- mle( minus_log_likelihood_geometric,
                        start = list(q = N / M),
                        method = "L-BFGS-B",
                        lower = c(0.000001),
                        upper = c(0.999999) )
  
  list(
    k = attributes(summary(mle_geometric))$coef[1],
    AIC = get_AIC(attributes(summary(mle_geometric))$m2logL, 2, N),
    distribution = function(x) {geometric_dist(x, attributes(summary(mle_geometric))$coef[1])})
}

######################### Altmann loglikelyhood ##################

loglikelyhood_altm <- function(degree_sequence){
  x <- degree_sequence$V1
  N <- length(x)
  M <- sum(x)
  C <- sum(sapply(x,function(j) sum(log(2:j))))
  
  minus_log_likelihood_altm <- function(p,q) {
    -sum(-N*log(sum((1:N)^(-p)*exp(-(1:N)*q))) + sum(log(x^(-p))) + sum(-x*q) ) 
  }
  
  mle_altman <- mle( minus_log_likelihood_altm,
                     start = list(p=2,q=0.5),
                     method = "L-BFGS-B",
                     lower = c(1.000001,0),
                     upper = c(2,1) )
  
  attributes(summary(mle_altman))$coef[2]
  
  list(
    p = attributes(summary(mle_altman))$coef[1],
    q = attributes(summary(mle_altman))$coef[2],
    AIC = get_AIC(attributes(summary(mle_altman))$m2logL, 1, N),
    distribution = function(x) {altm_dist(x, attributes(summary(mle_altman))$coef[1], attributes(summary(mle_altman))$coef[2])}  )
}

################## Main #########################

source = read.table("list.txt", 
                    header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
                    as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
)


barplot(degree_spectrum, main = "English (Bar Plot in normal scale)", xlab = "degree", ylab = "number of vertices", col=rgb(0.2,0.4,0.6,0.6) )
barplot(degree_spectrum, main = "English (Bar Plot in log-log scale)", xlab = "degree", ylab = "number of vertices", log = "xy", col=rgb(0.2,0.4,0.6,0.6))



# Without Altman
############################
cat("Results Without Altman\n")
for (x in 1:nrow(source)) {
  degree_sequence = read.table(source$file[x], header = FALSE)
  degree_spectrum = table(degree_sequence) / length(degree_sequence$V1)
  
  zeta_gamma2 <- loglikelyhood_zeta_gamma2(degree_sequence)
  zeta <- loglikelyhood_zeta(degree_sequence)
  zeta_truncated <- loglikelyhood_zeta_truncated(degree_sequence)
  poisson <- loglikelyhood_poisson(degree_sequence)
  geometric <- loglikelyhood_geometric(degree_sequence)
  
  
  AICs <- c(zeta_gamma2$AIC, zeta$AIC, zeta_truncated$AIC, poisson$AIC, geometric$AIC )
  AICs_diff <- AICs- min(AICs)
  
  
  
  rows = as.numeric(unlist(rownames(degree_spectrum)))
  
  cat("\n")
  cat("Data summary\n")
  write_summary(source$language[x], degree_sequence)
  
  cat("AIC diff\n")
  cat(source$language[x], " ")
  cat(round(AICs_diff))
  cat("\n")
  cat("Parameters\n")
  cat(source$language[x], " ")
  cat(zeta$gamma, zeta_truncated$gamma, zeta_truncated$kmax,poisson$k, geometric$k )
  cat("\n")
  cat("P(x>200)\n")
  cat(source$language[x],
      sum(degree_spectrum[rows>200]), 
      1 - sum(zeta_gamma2$distribution(seq(200))),
      1 - sum(zeta$distribution(seq(200))),
      1 - sum(zeta_truncated$distribution(seq(200))),
      1 - sum(poisson$distribution(seq(200))),
      1 - sum(geometric$distribution(seq(200))) ) 
  
  
  index = seq( max(degree_sequence$V1))
  
  n = log10(max(AICs_diff) + 1) + 1
  colors=hcl.colors(n + 3 , palette = "YlOrRd")[1:n]
  
  par(fig=c(0,0.85,0,1), mar=c(4,4,2,2))
  plot(rows,degree_spectrum, main = source$language[x], type = "p",
       xlab = "degree", ylab = "p", log = "xy", ylim=c(0.00001,1))
  aty <- axTicks(2)
  labels <- sapply(aty,function(i)
    as.expression(bquote( .(i)))
  )
  axis(2,at=aty,labels=labels)
  lines(zeta_gamma2$distribution(index), col = colors[log10(AICs_diff[1] + 1) + 1], lwd=2)
  lines(zeta$distribution(index), col = colors[log10(AICs_diff[2] + 1) + 1], lty=2, lwd=2)
  lines(zeta_truncated$distribution(index), col = colors[log10(AICs_diff[3] + 1) + 1], lty=3, lwd=2)
  lines(poisson$distribution(index), col = colors[log10(AICs_diff[4] + 1) + 1], lty=4, lwd=2)
  lines(geometric$distribution(index), col = colors[log10(AICs_diff[5] + 1) + 1], lty=5, lwd=2)
 
  
  legend("topright", 
         legend=c("Zeta (gamma=2)","Zeta", "Right Trunc. Zeta", "Poisson", "Geometric" ), 
         col=c(colors[log10(AICs_diff[1] + 1) + 1], colors[log10(AICs_diff[2] + 1) + 1], colors[log10(AICs_diff[3] + 1) + 1], colors[log10(AICs_diff[4] + 1) + 1], colors[log10(AICs_diff[5] + 1) + 1]  ), lty=1:6, cex=0.8, lwd=2)
  
  par(fig=c(0.85,1,0,1), new=TRUE, mar=c(1,1,2,3))
  image(y=(1:n),z=t(1:n), col=colors, axes=FALSE, main=expression(paste(Delta,"AIC")), cex.main=.8)
  aty <- axTicks(4)
  labels <- sapply(aty,function(i)
    as.expression(bquote( 10^.(i-1)))
  )
  axis(4, labels = labels, at=aty )
}


# With Altman
############################
cat("Results With Altman\n")
for (x in 1:nrow(source)) {
  degree_sequence = read.table(source$file[x], header = FALSE)
  degree_spectrum = table(degree_sequence) / length(degree_sequence$V1)
  
  zeta_gamma2 <- loglikelyhood_zeta_gamma2(degree_sequence)
  zeta <- loglikelyhood_zeta(degree_sequence)
  zeta_truncated <- loglikelyhood_zeta_truncated(degree_sequence)
  poisson <- loglikelyhood_poisson(degree_sequence)
  geometric <- loglikelyhood_geometric(degree_sequence)
  altm <- loglikelyhood_altm(degree_sequence)
  
  AICs <- c(zeta_gamma2$AIC, zeta$AIC, zeta_truncated$AIC, poisson$AIC, geometric$AIC, altm$AIC)
  AICs_diff <- AICs- min(AICs)
  
  
  
  rows = as.numeric(unlist(rownames(degree_spectrum)))
  
  cat("\n")
  cat("Data summary\n")
  write_summary(source$language[x], degree_sequence)
  
  cat("AIC diff\n")
  cat(source$language[x], " ")
  cat(round(AICs_diff))
  cat("\n")
  cat("Parameters\n")
  cat(source$language[x], " ")
  cat(zeta$gamma, zeta_truncated$gamma, zeta_truncated$kmax,poisson$k, geometric$k, altm$p, altm$q)
  cat("\n")
  cat("P(x>200)\n")
  cat(source$language[x],
      sum(degree_spectrum[rows>200]), 
      1 - sum(zeta_gamma2$distribution(seq(200))),
      1 - sum(zeta$distribution(seq(200))),
      1 - sum(zeta_truncated$distribution(seq(200))),
      1 - sum(poisson$distribution(seq(200))),
      1 - sum(geometric$distribution(seq(200))),
      1 - sum(altm$distribution(seq(200))))
  
  
  index = seq( max(degree_sequence$V1))
  
  n = log10(max(AICs_diff) + 1) + 1
  colors=hcl.colors(n + 3 , palette = "YlOrRd")[1:n]
  
  par(fig=c(0,0.85,0,1), mar=c(4,4,2,2))
  plot(rows,degree_spectrum, main = source$language[x], type = "p",
       xlab = "degree", ylab = "p", log = "xy", ylim=c(0.00001,1))
  aty <- axTicks(2)
  labels <- sapply(aty,function(i)
    as.expression(bquote( .(i)))
  )
  axis(2,at=aty,labels=labels)
  lines(zeta_gamma2$distribution(index), col = colors[log10(AICs_diff[1] + 1) + 1], lwd=2)
  lines(zeta$distribution(index), col = colors[log10(AICs_diff[2] + 1) + 1], lty=2, lwd=2)
  lines(zeta_truncated$distribution(index), col = colors[log10(AICs_diff[3] + 1) + 1], lty=3, lwd=2)
  lines(poisson$distribution(index), col = colors[log10(AICs_diff[4] + 1) + 1], lty=4, lwd=2)
  lines(geometric$distribution(index), col = colors[log10(AICs_diff[5] + 1) + 1], lty=5, lwd=2)
  lines(altm$distribution(index), col = colors[log10(AICs_diff[6] + 1) + 1], lty=6, lwd=2)
  
  legend("topright", 
         legend=c("Zeta (gamma=2)","Zeta", "Right Trunc. Zeta", "Poisson", "Geometric", "Altmann"), 
         col=c(colors[log10(AICs_diff[1] + 1) + 1], colors[log10(AICs_diff[2] + 1) + 1], colors[log10(AICs_diff[3] + 1) + 1], colors[log10(AICs_diff[4] + 1) + 1], colors[log10(AICs_diff[5] + 1) + 1], colors[log10(AICs_diff[6] + 1) + 1]), lty=1:6, cex=0.8, lwd=2)
  
  par(fig=c(0.85,1,0,1), new=TRUE, mar=c(1,1,2,3))
  image(y=(1:n),z=t(1:n), col=colors, axes=FALSE, main=expression(paste(Delta,"AIC")), cex.main=.8)
  aty <- axTicks(4)
  labels <- sapply(aty,function(i)
    as.expression(bquote( 10^.(i-1)))
  )
  axis(4, labels = labels, at=aty )
}

