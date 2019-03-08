library(haven)
library(Rcpp)

setwd("~/mywork/paper1/data")

set.seed(43)

sourceCpp("/home/bz/mywork/paper1/code/semibart/src/smmbart_arma.cpp")
source("/home/bz/mywork/paper1/code/semibart/R/semibart.R")

## should be same
mydata <- read_sas("bartdata.sas7bdat")
mydata2 <- read.csv("bartdata.csv")

## center fib4_ongoing_lag around 3.25
mydata2$fib4_ongoing_lag <- mydata2$fib4_ongoing_lag - 3.25
summary(mydata2$fib4_ongoing_lag)

## center other covariates around medians or meaningful values
median(mydata2$PC_AGE)  ## 50
mydata2$PC_AGE <- mydata2$PC_AGE - 50

median(mydata2$PC_YR_HAART) ## 2005
mydata2$PC_YR_HAART <- mydata2$PC_YR_HAART - 2005

median(mydata2$alt_ongoing_lag); summary(mydata2$alt_ongoing_lag) ## 50
mydata2$alt_ongoing_lag <- mydata2$alt_ongoing_lag - 50

median(mydata2$ast_ongoing_lag); summary(mydata2$ast_ongoing_lag) ## 50
mydata2$ast_ongoing_lag <- mydata2$ast_ongoing_lag - 50

median(mydata2$cd4_ongoing_lag); summary(mydata2$cd4_ongoing_lag) ## 200
mydata2$cd4_ongoing_lag <- mydata2$cd4_ongoing_lag - 200

median(mydata2$vl_log10_ongoing_lag); summary(mydata2$vl_log10_ongoing_lag) ## 4.75
mydata2$vl_log10_ongoing_lag <- mydata2$vl_log10_ongoing_lag - 4.75

x <- mydata2[ , c("tox_hep","PC_AGE","CO_ALC_ART","CO_DRUG_ART","PC_YR_HAART","alt_ongoing_lag",
                  "ast_ongoing_lag","cd4_ongoing_lag","dm_ongoing_lag","fib4_ongoing_lag",
                  "vl_log10_ongoing_lag","race_black","race_hisp","race_oth","bmilt25","bmi25_30")]

## x2 has no fib4_ongoing_lag
x2 <- mydata2[ , c("tox_hep","PC_AGE","CO_ALC_ART","CO_DRUG_ART","PC_YR_HAART","alt_ongoing_lag",
                   "ast_ongoing_lag","cd4_ongoing_lag","dm_ongoing_lag",
                   "vl_log10_ongoing_lag","race_black","race_hisp","race_oth","bmilt25","bmi25_30")]


## blips for 3 analyses
blip <- cbind(mydata2$tox_mit)
blip.em <- cbind(mydata2$tox_mit, mydata2$fib4_ongoing_lag, mydata2$tox_mit * mydata2$fib4_ongoing_lag)
blip.em2 <- cbind(mydata2$tox_mit, I(mydata2$fib4_ongoing_lag > 0), mydata2$tox_mit * I(mydata2$fib4_ongoing_lag > 0))

## 2 year binary indicator of death
y <- mydata2$death2yr


ntrees <- 50
nreps <- 20000
burnin <- 5000
power <- 2.0; base <- 0.95
sigdf <- 3; sigquant <- 0.90; k <- 2.0

smm_res <- semibart(x.train=as.matrix(x),
                    a.train=as.matrix(blip),
                    y.train=y,
                    sigest=NA,
                    sigdf=sigdf,
                    sigquant=sigquant,
                    k=k,
                    power=power,
                    base=base,
                    meanb=rep(0,dim(as.matrix(blip))[2]),
                    sigb=4,
                    ntree=ntrees,
                    ndpost=nreps,
                    numcut=rep(100,ncol(as.matrix(x))),
                    usequants=TRUE,
                    offset=0,
                    binarylink="probit",
                    verbose=TRUE,
                    printevery=100)

postscript("analysis1_tp_beta1.eps")
plot(1:nreps,smm_res$beta, type = "l", main = "Trace plot for mtNRTI coefficient", xlab = "Iteration", ylab = "Value")
dev.off()

mean(smm_res$beta[(burnin + 1):nreps])
quantile(smm_res$beta[(burnin + 1):nreps], probs = c(0.025, 0.5, 0.975))

probit1 <- glm(y ~ as.matrix(x) + as.matrix(blip), family = binomial("probit"))

smm_res2 <- semibart(x.train=as.matrix(x2),
                    a.train=as.matrix(blip.em),
                    y.train=y,
                    sigest=NA,
                    sigdf=sigdf,
                    sigquant=sigquant,
                    k=k,
                    power=power,
                    base=base,
                    meanb=rep(0,dim(as.matrix(blip.em))[2]),
                    sigb=4,
                    ntree=ntrees,
                    ndpost=nreps,
                    numcut=rep(100,ncol(as.matrix(x2))),
                    usequants=TRUE,
                    offset=0,
                    binarylink="probit",
                    verbose=TRUE,
                    printevery=250)



colMeans(smm_res2$beta[(burnin + 1):nreps, ])
quantile(smm_res2$beta[(burnin + 1):nreps, 1], probs = c(0.025, 0.5, 0.975))
quantile(smm_res2$beta[(burnin + 1):nreps, 3], probs = c(0.025, 0.5, 0.975))


probit2 <- glm(y ~ as.matrix(x2) + as.matrix(blip.em), family = binomial("probit"))


postscript("analysis2_tp_beta.eps")
par(mfrow=c(2,1)) 
plot(1:nreps,smm_res2$beta[ , 1], type='l', main="trace plot for mtNRTI coefficient", xlab = "Iteration", ylab = "Value")
plot(1:nreps,smm_res2$beta[ , 3], type='l', main="trace plot for effect modification coefficient (continuous)", xlab = "Iteration", ylab = "Value")
dev.off()


start <- Sys.time()
smm_res3 <- semibart(x.train=as.matrix(x2),
                     a.train=as.matrix(blip.em2),
                     y.train=y,
                     sigest=NA,
                     sigdf=sigdf,
                     sigquant=sigquant,
                     k=k,
                     power=power,
                     base=base,
                     meanb=rep(0,dim(as.matrix(blip.em2))[2]),
                     sigb=4,
                     ntree=ntrees,
                     ndpost=nreps,
                     numcut=rep(100,ncol(as.matrix(x2))),
                     usequants=TRUE,
                     offset=0,
                     binarylink="probit",
                     verbose=TRUE,
                     printevery=250)
end <- Sys.time()

probit3 <- glm(y ~ as.matrix(x2) + as.matrix(blip.em2), family = binomial("probit"))


colMeans(smm_res3$beta[(burnin + 1):nreps, ])
quantile(smm_res3$beta[(burnin + 1):nreps, 1], probs = c(0.025, 0.5, 0.975))
quantile(smm_res3$beta[(burnin + 1):nreps, 3], probs = c(0.025, 0.5, 0.975))

postscript("analysis3_tp_beta.eps")
par(mfrow=c(2,1)) 
plot(1:nreps,smm_res3$beta[ , 1], type='l', main = "trace plot for mtNRTI coefficient", xlab = "Iteration", ylab = "Value")
plot(1:nreps,smm_res3$beta[ , 3], type='l', main = "trace plot for effect modification coefficient (binary; > 3.25)", xlab = "Iteration", ylab = "Value")
dev.off()




#### summarize analyses with probit regression ####


summary(probit1)
confint(probit1)


summary(probit2)
confint(probit2)


summary(probit3)
confint(probit3)
