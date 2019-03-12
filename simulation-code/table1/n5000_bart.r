library(BayesTree)

load("../../data/n5000_sd1.Rdata")

truex1  <- 2
trueint <- -1
truex6  <- 2

bart.est <- numeric(500)
bart.in.ci <- numeric(500)
bart.ci.len <- numeric(500)


for(i in 1:500) {
  mydata <- n5000_sd1$dat4[i, ,]
  y      <- mydata[ ,  ncol(mydata)]
  x      <- mydata[ , -ncol(mydata)]
  nt     <- sum(x[,1])
  
  xa1 <- x[which(x[ ,1] == 1),]
  xa0 <- xa1 
  xa0[,1] <- 0 
  testX = rbind(xa1,xa0)
  
  mybart<-bart(y.train=y,x.train=x,x.test=testX,verbose=FALSE)
  tmp<-apply(mybart$yhat.test[,1:nt]-mybart$yhat.test[,(nt+1):(2*nt)],1,mean)
  tmp.ci <- quantile(tmp,probs=c(0.025,0.975))
  
  bart.est[i] <- mean(tmp)
  bart.in.ci[i] <- (tmp.ci[1] < truex1 & tmp.ci[2] > truex1)
  bart.ci.len[i] <- tmp.ci[2]-tmp.ci[1]
  cat("Done iteration: ",i,"\n")
}

save(bart.est,bart.in.ci,bart.ci.len,file="./output/BART_nl_n5000_sd1.Rdata")
