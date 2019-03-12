library(Rcpp)

sourceCpp("/home/zeldow/dissertation/smmbart/paper1/sims/finalsims/semibart/smmbart_arma.cpp")
source("/home/zeldow/dissertation/smmbart/paper1/sims/finalsims/semibart/semibart.R")

load("./data/n5000.Rdata")


truex1  <- 0.3
trueint <- -0.1
truex6  <- 0.1

ntot  <- 10000
nburn <- 2500

start <- 7 - 1
num <- 50

smm.res   <- matrix(NA, ncol = 3, nrow = num)
smm.in.ci <- matrix(NA, ncol = 3, nrow = num)

for( i in (num * start + 1):(num * start + num) ) {

  mydata <- n5000$dat6[i, ,]
  y      <- mydata[ ,  ncol(mydata)]
  x      <- mydata[ , -ncol(mydata)]

  set.seed(i)
  smm_res <- semibart(x.train = as.matrix(x[,-c(1, 2)]),
                      a.train = as.matrix(cbind(x[ ,1], x[ ,1]*x[ ,2], x[ ,2])),
                      y.train = y,
                      sigest = NA,
                      #sigest = 1,
                      sigdf  = 3,
                      sigquant = 0.90,
                      k = 2.0,
                      power = 2.0,
                      base = 0.95,
                      meanb = rep(0, 3),
                      sigb = 4,
                      ntree = 50,
                      ndpost = ntot,
                      numcut = 100,
                      usequants = TRUE,
                      offset = 0,
                      binarylink = "probit",
                      verbose = FALSE,
                      printevery = 1000)
smm.in.ci[i - num*start, 3] <- (quantile(smm_res$beta[(nburn + 1):ntot,3], probs = 0.025) < truex6 & quantile(smm_res$beta[(nburn + 1):ntot,3], probs = 0.975) > truex6)
smm.in.ci[i - num*start, 2] <- (quantile(smm_res$beta[(nburn + 1):ntot,2], probs = 0.025) < trueint & quantile(smm_res$beta[(nburn + 1):ntot,2], probs = 0.975) > trueint)

  smm.res[i - num*start,]   <- colMeans(smm_res$beta[(nburn + 1):ntot,])
  smm.in.ci[i - num*start,1] <- (quantile(smm_res$beta[(nburn + 1):ntot,1], probs = 0.025) < truex1 & quantile(smm_res$beta[(nburn + 1):ntot,1], probs = 0.975) > truex1)

  cat("DONE: ",i,"\n")
}

save(smm.res, smm.in.ci, file = "./output/n5000_bin7.Rdata")
