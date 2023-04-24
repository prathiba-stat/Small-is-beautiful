simulate_SCED <- function(T, beta, e.sigma, rho){
  t <- max(T)
  p <- nrow(beta)
  y <- matrix(NA, nrow(beta), t)
  e <- matrix(NA, 1, t)
  for (j in 1:p){
    e[1] <- rnorm(1, 0, e.sigma)
    for (i in 2:T[j]){
      e[i] <- rho*e[i - 1] + rnorm(1, 0, e.sigma)
    }
    y[j,1:T[j]] <- beta[j, 1] + beta[j, 2]*seq(0:(T[j] - 1)) + e[1:T[j]]
    trend <- beta[j, 1] + beta[j, 2]*((seq(1:T[j])) - 1)
  }
  y
}
T <- rep(8, 2)
beta <- matrix(c(2, 0, NA, NA), 2, 2, byrow = TRUE)
int.es <- 2
e.sigma <- 1
slope.es <- 1
beta[2, 1] <- int.es * e.sigma + beta[1, 1]
beta[2,2] <- (slope.es + (beta[1, 1] + (T[1] + (T[2]/2))  * beta[1,2]) - 
                beta[2, 1])/(T[1] + (T[2]/2)) 
y <- simulate_SCED(T, 
                   beta, 
                   e.sigma, .5)
plot_SSD_Sim(y, "simulate.jpg")
