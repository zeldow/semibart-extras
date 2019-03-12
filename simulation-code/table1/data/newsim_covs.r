library(mvnfast)

p1 <- 20 # no. of continuous covariates
p2 <- 5 # no. of binary

makedata <- function(nsim, n, p1, p2, sd) {

  dat1 <- array(0, dim = c(nsim, n, p1 + p2 + 1) )
  dat2 <- array(0, dim = c(nsim, n, p1 + p2 + 1) )
  dat3 <- array(0, dim = c(nsim, n, p1 + p2 + 1) )
  dat4 <- array(0, dim = c(nsim, n, p1 + p2 + 1) )
  dat5 <- array(0, dim = c(nsim, n, p1 + p2 + 1) )
  dat6 <- array(0, dim = c(nsim, n, p1 + p2 + 1) )
  dat7 <- array(0, dim = c(nsim, n, p1 + p2 + 1) )
  dat8 <- array(0, dim = c(nsim, n, p1 + p2 + 1) )
  dat9 <- array(0, dim = c(nsim, n, p1 + p2 + 1) )
  dat10 <- array(0, dim = c(nsim, n, p1 + p2 + 1) )


  for (k in 1:nsim) {
    cov <- matrix(0, p1, p1)
    diag(cov) <- rep(1, p1)

    for(i in 1:5) {
      for(j in 1:5) {
        if (i != j) cov[i,j] <- cov[j, i] <- 0.20
      }
    }

    for(i in 6:10) {
      for(j in 6:10) {
        if (i != j) cov[i,j] <- cov[j, i] <- 0.15
      }
    }

    for(i in 11:15) {
      for(j in 11:15) {
        if (i != j) cov[i,j] <- cov[j, i] <- 0.10
      }
    }

    for(i in 16:20) {
      for(j in 16:20) {
        if (i != j) cov[i,j] <- cov[j, i] <- 0.05
      }
    }

    mu <- c(rep(2.0, 5), rep(1.5, 5), rep(1.0, 5), rep(0.0, 5))

    cont.covs <- rmvn(n, mu, sigma = cov)

    bin.covs <- cbind(rbinom(n, size = 1, prob = 0.25), 
                      rbinom(n, size = 1, prob = 0.5), 
                      rbinom(n, size = 1, prob = 0.5),
                      rbinom(n, size = 1, prob = 0.75),
                      rbinom(n, size = 1, prob = 0.75))
    
    x <- cbind(bin.covs, cont.covs)

    ## linear binary - single
    prob.y1 <- pnorm(0.1 + 0.3 * x[ , 1]  + 0.04 * x[ , 6] - 0.02 * x[ , 7] + 0.1 * x[ , 2] - 0.03 * x[ , 10] + 0.04 * x[ , 9])
    y1      <- rbinom(n, size = 1, prob = prob.y1)


    ## nonlinear binary - single
    prob.y2 <- pnorm(0.1 + 0.3 * x[ , 1]  + 0.04 * x[ , 6] - sin( pi / 4 * x[ , 2] * x[ , 7] ) + exp(x[,7] / 10) -
                     0.02 * x[, 2] * x[ , 9] * x[ ,10])
    y2      <- rbinom(n, size = 1, prob = prob.y2)
  
    
    ## linear continuous - single
    mu3 <- 1 + 2 * x[ , 1] + 2 * x[ , 6] - 0.5 * x[ , 7] + 2 * x[ , 5] - 1.5 * x[ , 10] - 0.5 * x[ , 8]
    y3  <- rnorm(n, mu3 , sd = sd)


    ## nonlinear continuous  single
    mu4 <- 1 + 2 * x[ , 1] + 2 * x[ , 6] + sin( pi * x[ , 2] * x[ ,7] ) - 2 * exp( x[ , 3] * x[ , 5] ) +
           log( abs( cos ( pi / 2 * x[ , 8] ) ) ) - 1.8 * cos( x[ , 9]) + 3 * x[ , 3] * abs(x[ ,7]) ^ 1.5
    y4  <- rnorm(n, mu4 , sd = sd)

    ## linear binary - multi
    prob.y5 <- pnorm(0.1 + 0.3 * x[ , 1] - 0.1 * x[ , 1] * x[ , 6] + 0.04 * x[ , 6] - 0.02 * x[ , 7] + 0.1 * x[ , 2] - 0.03 * x[ , 10] + 0.04 * x[ , 9])
    y5      <- rbinom(n, size = 1, prob = prob.y5)


    ## nonlinear binary - multi
    prob.y6 <- pnorm(0.1 + 0.3 * x[ , 1] - 0.1 * x[ , 1] * x[ , 6] + 0.04 * x[ , 6] - sin( pi / 4 * x[ , 2] * x[ , 7] ) + exp(x[,7] / 10) -
                     0.02 * x[, 2] * x[ , 9] * x[ ,10])
    y6      <- rbinom(n, size = 1, prob = prob.y6)
  
    
    ## linear continuous - multi
    mu7 <- 1 + 2 * x[ , 1] - 1 * x[ , 1] * x[ , 6] + 2 * x[ , 6] - 0.5 * x[ , 7] + 2 * x[ , 5] - 1.5 * x[ , 10] - 0.5 * x[ , 8]
    y7  <- rnorm(n, mu7 , sd = sd)


    ## nonlinear continuous -multi
    mu8 <- 1 + 2 * x[ , 1] - 1 * x[ , 1] * x[ , 6] + 2 * x[ , 6] + sin( pi * x[ , 2] * x[ ,7] ) - 2 * exp( x[ , 3] * x[ , 5] ) +
           log( abs( cos ( pi / 2 * x[ , 8] ) ) ) - 1.8 * cos( x[ , 9]) + 3 * x[ , 3] * abs(x[ ,7]) ^ 1.5
    y8  <- rnorm(n, mu8 , sd = sd)
  
    ## nonlinear binary - multi mispec
    prob.y9 <- pnorm(0.1 + 0.3 * x[ , 1] - 0.1 * x[ , 1] * x[ , 6] + sin( 0.04 * x[ , 6] ) - sin( pi / 4 * x[ , 2] * x[ , 7] ) + exp(x[,7] / 10) -
                     0.02 * x[, 2] * x[ , 9] * x[ ,10])
    y9      <- rbinom(n, size = 1, prob = prob.y9)
    
    ## nonlinear continuous -multi mispec
    mu10 <- 1 + 2 * x[ , 1] - 1 * x[ , 1] * x[ , 6] + sin(  2 * x[ , 6] ) + sin( pi * x[ , 2] * x[ ,7] ) - 2 * exp( x[ , 3] * x[ , 5] ) +
           log( abs( cos ( pi / 2 * x[ , 8] ) ) ) - 1.8 * cos( x[ , 9]) + 3 * x[ , 3] * abs(x[ ,7]) ^ 1.5
    y10  <- rnorm(n, mu10 , sd = sd)
    
    mydat1      <- cbind(x, y1)
    dat1[k, , ] <- mydat1

    mydat2      <- cbind(x, y2)
    dat2[k, , ] <- mydat2
    
    mydat3      <- cbind(x, y3)
    dat3[k, , ] <- mydat3
    
    mydat4      <- cbind(x, y4)
    dat4[k, , ] <- mydat4

    mydat5      <- cbind(x, y5)
    dat5[k, , ] <- mydat5

    mydat6      <- cbind(x, y6)
    dat6[k, , ] <- mydat6
    
    mydat7      <- cbind(x, y7)
    dat7[k, , ] <- mydat7
    
    mydat8      <- cbind(x, y8)
    dat8[k, , ] <- mydat8

    mydat9      <- cbind(x, y9)
    dat9[k, , ] <- mydat9
    
    mydat10      <- cbind(x, y10)
    dat10[k, , ] <- mydat10
  }

  return(list( dat1 = dat1, dat2 = dat2, dat3 = dat3, dat4 = dat4, dat5 = dat5, 
               dat6 = dat6, dat7 = dat7, dat8 = dat8, dat9 = dat9, dat10 = dat10) )

}

set.seed(44)

n250_sd01  <- makedata(500, 250, 20, 5, 0.1)
n1000_sd01 <- makedata(500, 1000, 20, 5, 0.1)
n5000_sd01 <- makedata(500, 5000, 20, 5, 0.1)

n250_sd1  <- makedata(500, 250, 20, 5, 1)
n1000_sd1 <- makedata(500, 1000, 20, 5, 1)
n5000_sd1 <- makedata(500, 5000, 20, 5, 1)

n250_sd2  <- makedata(500, 250, 20, 5, 2)
n1000_sd2 <- makedata(500, 1000, 20, 5, 2)
n5000_sd2 <- makedata(500, 5000, 20, 5, 2)

n250_sd3  <- makedata(500, 250, 20, 5, 3)
n1000_sd3 <- makedata(500, 1000, 20, 5, 3)
n5000_sd3 <- makedata(500, 5000, 20, 5, 3)
n500_sd1 <- makedata(500, 500, 20, 5, 1)

save(n500_sd1, file = "n500_sd1.Rdata")

# save(n250_sd01, file = "n250_sd01.Rdata")
# save(n250_sd1, file = "n250_sd1.Rdata")
# save(n250_sd2, file = "n250_sd2.Rdata")
# save(n250_sd3, file = "n250_sd3.Rdata")
# 
# save(n1000_sd01, file = "n1000_sd01.Rdata")
# save(n1000_sd1, file = "n1000_sd1.Rdata")
# save(n1000_sd2, file = "n1000_sd2.Rdata")
# save(n1000_sd3, file = "n1000_sd3.Rdata")
# 
# save(n5000_sd01, file = "n5000_sd01.Rdata")
# save(n5000_sd1, file = "n5000_sd1.Rdata")
# save(n5000_sd2, file = "n5000_sd2.Rdata")
# save(n5000_sd3, file = "n5000_sd3.Rdata")
