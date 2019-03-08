Purpose
-------

This folder houses some ancillary simulations we did to assess bias in
models that are part-tree, part-linear, such as semi-BART. These
simulations were done to help with our understanding and are not
comprehensive or conclusive.

Setup
-----

We simulate one covariate *x* which is iid standard normal.
Independently, we simulate a binary treatment *a* from iid
Bernoulli(0.5). The outcome is generated with a normal error term with
mean in the form of
*μ*(*a*, *x*) = *ψ*<sub>1</sub>*a* + *ψ*<sub>2</sub>*a**x* + ϖ(*x*)
where ϖ(*x*) is a regression tree dependent only on *x*.

These simulations use three structures for ϖ(*x*):

-   Very basic tree (one split point): tree\_simplest.R
-   Basic tree (two split points): tree\_simple.R
-   Complicated tree (many split points): tree\_complex.R

Estimation
----------

Estimation is ad hoc and not intended to be representative of a real
data analysis. We use a backfitting algorithm to estimate the tree part
and the linear part. The model is correctly specified. Basically we
estimate a tree and get estimates. We get residuals by subtracting off
the fit of the tree, then we estimate the linear part. We then subtract
off the linear fit and estimate the tree again, continuing this process
until convergence.
