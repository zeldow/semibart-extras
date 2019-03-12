library(mvnfast)

expit <- function(x) exp(x) / (1 + exp(x))

makedata <- function(nsim, n) {

  dat6 <- array(0, dim = c(nsim, n, 30 + 2) )
  dat8 <- array(0, dim = c(nsim, n, 30 + 2) )


  for (k in 1:nsim) {
    p <- 30
    mu <- rep(0, p)
    sig <- matrix(0, nrow = p, ncol = p)
    diag(sig) <- 1
    rho <- 0.5
    for(i in 1:(p-1)) {
      for(j in (i+1):p) {
        kk <- j - i
        sig[i, j] <- sig[j, i] <- rho^kk
      }
    }

    x <- rmvn(n, mu, sigma = sig)

   #  bin.covs <- cbind(rbinom(n, size = 1, prob = 0.25), 
   #                    rbinom(n, size = 1, prob = 0.5), 
   #                    rbinom(n, size = 1, prob = 0.5),
   #                    rbinom(n, size = 1, prob = 0.75),
   #                    rbinom(n, size = 1, prob = 0.75))
    
    prob.a <- expit(0.1 + 0.2 * x[,1] - sin(x[,3])/3 - 0.1 * x[,22])
    a <- rbinom(n, 1, prob.a)
    x <- cbind(a, x)

#     ## linear binary - single
#     prob.y1 <- pnorm(0.1 + 0.3 * x[ , 1]  + 0.04 * x[ , 6] - 0.02 * x[ , 7] + 0.1 * x[ , 2] - 0.03 * x[ , 10] + 0.04 * x[ , 9])
#     y1      <- rbinom(n, size = 1, prob = prob.y1)
# 
# 
#     ## nonlinear binary - single
#     prob.y2 <- pnorm(0.1 + 0.3 * x[ , 1]  + 0.04 * x[ , 6] - sin( pi / 4 * x[ , 2] * x[ , 7] ) + exp(x[,7] / 10) -
#                      0.02 * x[, 2] * x[ , 9] * x[ ,10])
#     y2      <- rbinom(n, size = 1, prob = prob.y2)
#   
#     
#     ## linear continuous - single
#     mu3 <- 1 + 2 * x[ , 1] + 2 * x[ , 6] - 0.5 * x[ , 7] + 2 * x[ , 5] - 1.5 * x[ , 10] - 0.5 * x[ , 8]
#     y3  <- rnorm(n, mu3 , sd = sd)
# 
# 
#     ## nonlinear continuous  single
#     mu4 <- 1 + 2 * x[ , 1] + 2 * x[ , 6] + sin( pi * x[ , 2] * x[ ,7] ) - 2 * exp( x[ , 3] * x[ , 5] ) +
#            log( abs( cos ( pi / 2 * x[ , 8] ) ) ) - 1.8 * cos( x[ , 9]) + 3 * x[ , 3] * abs(x[ ,7]) ^ 1.5
#     y4  <- rnorm(n, mu4 , sd = sd)
# 
#     ## linear binary - multi
#     prob.y5 <- pnorm(0.1 + 0.3 * x[ , 1] - 0.1 * x[ , 1] * x[ , 6] + 0.04 * x[ , 6] - 0.02 * x[ , 7] + 0.1 * x[ , 2] - 0.03 * x[ , 10] + 0.04 * x[ , 9])
#     y5      <- rbinom(n, size = 1, prob = prob.y5)
# 
# 
     ## nonlinear binary - multi
     prob.y6 <- pnorm(0.1 + 0.3 * x[ , 1] - 0.1 * x[ , 1] * x[ , 2] + 0.1 * x[ , 2] - sin( pi / 4 * x[ , 22] * x[ , 7] ) + exp(x[,7] / 5)*x[,11]/4 -
                      0.12 * x[, 22] * x[ , 9] * x[ ,10] + 0.05*x[,8] * x[,10]*x[,11]^2)
     y6      <- rbinom(n, size = 1, prob = prob.y6)
#   
#     
#     ## linear continuous - multi
#     mu7 <- 1 + 2 * x[ , 1] - 1 * x[ , 1] * x[ , 6] + 2 * x[ , 6] - 0.5 * x[ , 7] + 2 * x[ , 5] - 1.5 * x[ , 10] - 0.5 * x[ , 8]
#     y7  <- rnorm(n, mu7 , sd = sd)
# 

    ## nonlinear continuous -multi
    mu8 <- 1 + 2 * x[ , 1] - 1 * x[ , 1] * x[ , 2] + 2 * x[ , 2] + sin( pi * x[ , 22] * x[ ,7] ) - 1 * exp( x[ , 6]/5 * x[ , 5] ) +
           log( abs( cos ( pi / 2 * x[ , 8] ) ) ) - 1.8 * cos( x[ , 9]) + 0.2 * x[ , 11] * abs(x[ ,7]) ^ 1.5
    y8  <- rnorm(n, mu8 , sd = 1)
  
#     ## nonlinear binary - multi mispec
#     prob.y9 <- pnorm(0.1 + 0.3 * x[ , 1] - 0.1 * x[ , 1] * x[ , 6] + sin( 0.04 * x[ , 6] ) - sin( pi / 4 * x[ , 2] * x[ , 7] ) + exp(x[,7] / 10) -
#                      0.02 * x[, 2] * x[ , 9] * x[ ,10])
#     y9      <- rbinom(n, size = 1, prob = prob.y9)
#     
#     ## nonlinear continuous -multi mispec
#     mu10 <- 1 + 2 * x[ , 1] - 1 * x[ , 1] * x[ , 6] + sin(  2 * x[ , 6] ) + sin( pi * x[ , 2] * x[ ,7] ) - 2 * exp( x[ , 3] * x[ , 5] ) +
#            log( abs( cos ( pi / 2 * x[ , 8] ) ) ) - 1.8 * cos( x[ , 9]) + 3 * x[ , 3] * abs(x[ ,7]) ^ 1.5
#     y10  <- rnorm(n, mu10 , sd = sd)
#     
#     mydat1      <- cbind(x, y1)
#     dat1[k, , ] <- mydat1
# 
#     mydat2      <- cbind(x, y2)
#     dat2[k, , ] <- mydat2
#     
#     mydat3      <- cbind(x, y3)
#     dat3[k, , ] <- mydat3
#     
#     mydat4      <- cbind(x, y4)
#     dat4[k, , ] <- mydat4
# 
#     mydat5      <- cbind(x, y5)
#     dat5[k, , ] <- mydat5
# 
    mydat6      <- cbind(x, y6)
    dat6[k, , ] <- mydat6
#     
#     mydat7      <- cbind(x, y7)
#     dat7[k, , ] <- mydat7
    
    mydat8      <- cbind(x, y8)
    dat8[k, , ] <- mydat8

#     mydat9      <- cbind(x, y9)
#     dat9[k, , ] <- mydat9
#     
#     mydat10      <- cbind(x, y10)
#     dat10[k, , ] <- mydat10

  }

  return(list( dat6 = dat6, dat8 = dat8) )

}

set.seed(44)


n250  <- makedata(500, 250)
n500  <- makedata(500, 500)
n1000 <- makedata(500, 1000)
n5000 <- makedata(500, 5000)


save(n250, file = "n250.Rdata")
save(n500, file = "n500.Rdata")
save(n1000, file = "n1000.Rdata")
save(n5000, file = "n5000.Rdata")
