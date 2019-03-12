library(mgcv)

load("./data/n500.Rdata")


truex1  <- 2
trueint <- -1
truex6  <- 2


gam.res   <- matrix(NA, ncol = 3, nrow = 500)
gam.in.ci <- matrix(NA, ncol = 3, nrow = 500)


for( i in 1:500) {

  mydata <- n500$dat8[i, ,]
  y      <- mydata[ ,  ncol(mydata)]
  x      <- mydata[ , -ncol(mydata)]

  g <- gam(y ~ x[,1] + x[,2] + x[,1]*x[,2] + s(x[,3]) + s(x[,4]) +
               s(x[,5]) + s(x[,6]) + s(x[,7]) + s(x[,8]) + s(x[,9]) +
               s(x[,10]) + s(x[,11]) + s(x[,12]) + s(x[,13]) + s(x[,14]) +
               s(x[,15]) + s(x[,16]) + s(x[,17]) + s(x[,18]) + s(x[,19]) +
               s(x[,20]) + s(x[,21]) + s(x[,22]) + s(x[,23]) + s(x[,24]) +
               s(x[,25]) + s(x[,26]) + s(x[,27]) + s(x[,28]) + s(x[,29]) + 
               s(x[,30]) + s(x[,31]))

  gam.res[i, ] <- g$coefficients[c(2,4,3)]
  gam.in.ci[i, 1] <- (g$coefficients[2] - qnorm(0.975)*summary(g)$se[2] < truex1 & 
                      g$coefficients[2] + qnorm(0.975)*summary(g)$se[2] > truex1)

  gam.in.ci[i, 2] <- (g$coefficients[4] - qnorm(0.975)*summary(g)$se[4] < trueint & 
                      g$coefficients[4] + qnorm(0.975)*summary(g)$se[4] > trueint)
  
  gam.in.ci[i, 3] <- (g$coefficients[3] - qnorm(0.975)*summary(g)$se[3] < truex6 & 
                      g$coefficients[3] + qnorm(0.975)*summary(g)$se[3] > truex6)

  cat("Iteration: ", i, " of 500\n")
}

save(gam.res, gam.in.ci, file = "./output/n500_cont_gam.Rdata")
