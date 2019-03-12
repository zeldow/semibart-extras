#library(Rcpp)

#sourceCpp("/home/zeldow/dissertation/smmbart/paper1/sims/finalsims/code/smmbart_arma.cpp")
#source("/home/zeldow/dissertation/smmbart/paper1/sims/finalsims/code/semibart.R")

load("../../../data/n5000_sd1.Rdata")

truex1  <- 2
trueint <- -1
truex6  <- 2

lm.res   <- matrix(NA, ncol = 1, nrow = 500)
lm.in.ci <- matrix(NA, ncol = 1, nrow = 500)

for( i in 1:500) {

  mydata <- n5000_sd1$dat4[i, ,]
  y      <- mydata[ ,  ncol(mydata)]
  x      <- mydata[ , -ncol(mydata)]

  reg    <- lm(y ~ x)
  lm.res[i,]   <- reg$coefficients[c("x1")]
  lm.in.ci[i,1] <- (confint(reg, "x1")[1] < truex1 & confint(reg, "x1")[2] > truex1)

}

save(lm.res, lm.in.ci, file = "./output/n5000_cont_nl_sd1.Rdata")
