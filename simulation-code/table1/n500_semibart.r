library(Rcpp)

sourceCpp("/home/zeldow/dissertation/smmbart/paper1/sims/finalsims/semibart/smmbart_arma.cpp")
source("/home/zeldow/dissertation/smmbart/paper1/sims/finalsims/semibart/semibart.R")

load("./data/n500_sd1.Rdata")

truex1  <- 2
trueint <- -1
truex6  <- 2

ntot  <- 10000
nburn <- 2500

start <- 1 - 1
num <- 25

smm.res   <- matrix(NA, ncol = 1, nrow = num)
smm.in.ci <- matrix(NA, ncol = 1, nrow = num)

for( i in (num * start + 1):(num * start + num) ) {

  mydata <- n500_sd1$dat4[i, ,]
  y      <- mydata[ ,  ncol(mydata)]
  x      <- mydata[ , -ncol(mydata)]

  set.seed(i)
  smm_res <- semibart(x.train = as.matrix(x[,-1]),
                      a.train = as.matrix(x[ ,1]),
                      y.train = y,
                      sigest = NA,
                      sigdf  = 3,
                      sigquant = 0.90,
                      k = 2.0,
                      power = 2.0,
                      base = 0.95,
                      meanb = 0,
                      sigb = 4,
                      ntree = 50,
                      ndpost = ntot,
                      numcut = 100,
                      usequants = TRUE,
                      offset = -9999.0,
                      binarylink = "probit",
                      verbose = TRUE,
                      printevery = 1000)

  smm.res[i - num*start,]   <- mean(smm_res$beta[(nburn + 1):ntot])
  smm.in.ci[i - num*start,1] <- (quantile(smm_res$beta[(nburn + 1):ntot], probs = 0.025) < truex1 & quantile(smm_res$beta[(nburn + 1):ntot], probs = 0.975) > truex1)

  cat("DONE: ",i,"\n")
}

save(smm.res, smm.in.ci, file = "./output/n500_sb.Rdata")
