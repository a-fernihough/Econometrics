# Estimation of Treatment Effects Without an Exclusion Restriction

The script `bmate.R` estimates the bias-minimizing treatment effects as in Millimet and Tchernis (2013). It will thus produce results that are compatible with the Stata command `bmte`. Please read the Millimet and Tchernis (2013) and the McCarthy, Millimet and Tchernis (2014) for further details.
***
### References

Millimet, D.L. and Tchernis, R. (2013). Estimation of Treatment Effects Without an Exclusion Restriction: with an Application to the Analysis of the School Breakfast Program. *Journal of Applied Econometrics*, 28: 982-1017. https://doi.org/10.1002/jae.2286

McCarthy, I. & Millimet, D.L. and Tchernis, R. (2014). The Bmte Command: Methods for the Estimation of Treatment Effects when Exclusion Restrictions are Unavailable. The Stata Journal: Promoting communications on statistics and Stata. 14. 670-683. https://doi.org/10.1177/1536867X1401400311 

Klein, R. and Vella, F. (2009), A Semiparametric Model for Binary Response and Continuous Outcomes Under Index Heteroscedasticity. *Journal of Applied Econometrics*, 24(5), 735-762. https://doi.org/10.1002/jae.1064

***
### Example Code
```
library(MASS) # for multivariate normals
library(glmx) # for hetero probits
require(AER) # for ivreg function
library(devtools) # to load raw code

# bmte functions
source_url("https://raw.githubusercontent.com/a-fernihough/Econometrics/main/BMTE/bmate.R")

# simulate data
set.seed(112)
n = 10000
x1 = runif(n,-1,1) 
x2 = runif(n,-1,1)
hx = 0.5*(x1-x2) + 0.5*(x1^2-x2^2) + 2*x1*x2
Sigma = matrix(c(1,1,-0.25,1,1,-0.25,-0.25,-0.25,1),3,3)
j1 = mvrnorm(n,rep(0,3),Sigma)
t = ifelse(0.5 + hx - j1[,3] >0, 1, 0) 
y = ifelse(t==1,1+hx+j1[,1],hx+j1[,2])
data = data.frame(cbind(y,x1,x2,t))  

# get bmate with defaults
bmate(y ~ t | x1 + x2 + I(x1^2) + I(x2^2) + I(x1*x2), 
      data = data, theta=0.25, digits= 3)

# now explore estimates with Klein-Vella method
set.seed(611165)
n = 10000
x1 = runif(n,-1,1) 
x2 = runif(n,-1,1)
hx = 0.5*(x1-x2) + 0.5*(x1^2-x2^2) + 2*x1*x2
Sigma = matrix(c(1,1,-0.25,1,1,-0.25,-0.25,-0.25,1),3,3)
j1 = mvrnorm(n,rep(0,3),Sigma)
j1[,3] = (1+0.45*(x1+x2))*j1[,3]
t = ifelse(0.5 + hx - j1[,3] >0, 1, 0) 
y = ifelse(t==1,1+hx+j1[,1],hx+j1[,2])
data = data.frame(cbind(y,x1,x2,t))  

# Klein-Vella
kvate(y ~ t | x1 + x2 + I(x1^2) + I(x2^2) + I(x1*x2), 
      data=data)
```
