library(mvnfast)
library(rpart)

### regression simulation with tree
nsims <- 500
ests <- matrix(NA, nrow = nsims, ncol = 2)
for(j in 1:nsims) {
n <- 5000
x <- rnorm(n)
a <- rbinom(n,1,0.5)
mp <- 2*a + 1*a*x
y <- numeric(n)

sd <- 1


### set up outcome as a tree (simple version)
for(i in 1:n) {
  if(x[i] < -1) y[i] <- rnorm(1, mean = mp[i] + 2, sd = sd)
  else if(x[i] > 1) y[i] <- rnorm(1, mean = mp[i] + 4, sd = sd)
  else y[i] <- rnorm(1, mean = mp[i] + 3, sd = sd)
}
X <- cbind(a, a*x)

p <- ncol(X)


### Fit Bayesian model
# priors
alp <- 1
bet <- 1

Sigma <- 100 * diag(p)
Sig.vec <- diag(Sigma)
beta0 <- rep(0, p)


# posterior
alp.post <- alp + n / 2


# initial values
beta <- rep(0, p)
sig <- 1



ngibbs <- 500

beta.reps <- matrix(NA, nrow=ngibbs, ncol=p)
sig.reps <- numeric(ngibbs)


## all together
for(i in 1:ngibbs) {
  
  ## do fit of tree
  ystar <- y - X%*%beta  # subtract off fit of regression
  tree <- rpart(ystar ~ x)
  
  
  ## regression
  ystar <- y - predict(tree)
  
  v<-solve(solve(Sigma)+(1/sig)*t(X)%*%X)
  m<-v%*%((1/Sig.vec)*beta0+(1/sig)*t(X)%*%ystar)  # more generally, diag(t0mat)%*%m0
  beta.reps[i,]<- beta <-c(rmvn(1,m,v))
  bet.post <- bet + sum((ystar-X%*%beta)^2)/2
  
  
  sig.reps[i] <- sig <- 1/rgamma(1,alp.post,bet.post)
  
}

ests[j,] <- colMeans(beta.reps)
}

save(ests, file = "tree_simple.Rdata")
