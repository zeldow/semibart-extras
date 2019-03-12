
load("./data/n5000.Rdata")

truex1  <- 2
trueint <- -1
truex6  <- 2

lm.res   <- matrix(NA, ncol = 3, nrow = 500)
lm.in.ci <- matrix(NA, ncol = 3, nrow = 500)

for( i in 1:500) {

  mydata <- n5000$dat8[i, ,]
  y      <- mydata[ ,  ncol(mydata)]
  x      <- mydata[ , -ncol(mydata)]

  reg    <- lm(y ~ x + x[,1]:x[,2])
  lm.res[i,]   <- reg$coefficients[c("x1", "x[, 1]:x[, 2]", "x2")]
  lm.in.ci[i,1] <- (confint(reg, "x1")[1] < truex1 & confint(reg, "x1")[2] > truex1)
  lm.in.ci[i,2] <- (confint(reg, "x[, 1]:x[, 2]")[1] < trueint & confint(reg, "x[, 1]:x[, 2]")[2] > trueint)
  lm.in.ci[i,3] <- (confint(reg, "x2")[1] < truex6 & confint(reg, "x2")[2] > truex6)

}

save(lm.res, lm.in.ci, file = "./output/n5000_cont_reg.Rdata")
