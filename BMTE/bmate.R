#rm(list=ls())

library(glmx)
require(AER)

kvest = function(formula, data, digits=3, maxit=500, method= "BFGS"){
  # formula parsing
  outcome = as.character(formula)[2]
  formula = as.formula(gsub("\\|","~",as.character(formula)[3]))
  
  # run hetprob model
  h1 = try(hetglm(formula, family = binomial(link = "probit"), data, maxit=maxit, method=method))                  
  if(class(h1) == "try-error" | h1$converged==0){
    stop("Heteroskedastic probit model did not converge. Increasing maxit or using an alternative method may help")
  }  
  pf = paste(as.character(formula[2]), unlist(strsplit(as.character(formula[3]),"\\|"))[1], sep=" ~ ")
  p1 = glm(pf, family = binomial(link = "probit"), data)
  #save lr test
  lrhet = suppressWarnings(lrtest(p1, h1))
  # get predictions
  data$phstar = predict(h1, newdata=data)
    
  # formula parsing for non-heterogenous effects
  f1 = paste(as.character(formula[2]), unlist(strsplit(as.character(formula[3]),"\\|"))[1], sep=" + ")
  f2 = paste("phstar", unlist(strsplit(as.character(formula[3]),"\\|"))[1], sep=" + ")
  f3 = formula(paste(outcome, paste(f1 , f2,sep=" | "),sep= " ~ "))
  iv1 = ivreg(f3, data=data)
  ate1 = coef(iv1)[names(coef(iv1))==as.character(formula[2])]
    
  # formula parsing for heterogenous effects
  v1 = gsub(" ", "", unlist(strsplit(f1,"\\+")))
  v2 = paste(v1[1],v1[-1],sep=":")
  f4 = paste(f1, paste(v2,collapse=" + "),sep=" + ")
    
  v3 = gsub(" ", "", unlist(strsplit(f2,"\\+")))
  v4 = paste(v3[1],v3[-1],sep=":")
  f5 = paste(f2, paste(v4,collapse=" + "),sep=" + ")
  f6 = formula(paste(outcome, paste(f4 , f5,sep=" | "),sep= " ~ "))
  iv2 = ivreg(f6, data=data)
  
  # ate = ate + treatment interactions
  ate2a = coef(iv2)[names(coef(iv2))==as.character(formula[2])]
  # get a model matrix so the factors are as numerics
  mat1 = model.matrix(formula(paste("~",paste(v1[-1], collapse="+"),sep="")), data)
  mat1 = mat1[,colnames(mat1)!=c("(Intercept)")]
  # fill in NA coefs with zeros
  coef1 = iv2$coef[paste(v1[1],colnames(mat1),sep=":")]
  coef1[is.na(coef1)] = 0 
  ate2b = colMeans(mat1) %*% coef1
  ate2 = ate2a + ate2b
    
  # gather results
  ate = as.numeric(sprintf(paste("%.",paste(digits,"f",sep=""),sep=""),c(ate1,ate2)))
  results = list(cate=ate[1], hate=ate[2], thetprob=h1, lrhet=lrhet)  
  return(results)
}

kvate = function(x, ...) UseMethod("kvate")

kvate.default = function(formula, data, digits=3, maxit=500, method= "BFGS"){
  est = kvest(formula, data, digits, maxit, method)
  est$call = match.call() 
  class(est) = "kvate"
  est
}

print.kvate = function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nConstant Average Treatment Effect    :", x$cate)
  cat("\nHeterogenous Average Treatment Effect:",x$hate)
  cat("\n\nTest for Heteroscedasticity:\n\n")
  print(x$lrhet)
}

bmest = function(formula, data, theta=0.25, digits=3){
  
  if(theta>1 | theta<0){stop("Please choose a value of theta between 0 and 1")}
  
  # parse formula
  outcome = as.character(formula)[2]
  formula = as.formula(gsub("\\|","~",as.character(formula)[3]))  
  f1 = paste(as.character(formula)[2],as.character(formula)[3],sep=" + ")
  f2 = formula(paste(outcome, f1, sep=" ~ "))
  
  # ols non-heterogenous
  ols = lm(f2,data)
  ateols = coef(ols)[names(coef(ols))==as.character(formula)[2]]
  
  # ols hetero formula parsing
  v1 = gsub(" ", "",unlist(strsplit(as.character(formula)[3],"\\+")))
  f3 = paste(paste(as.character(formula)[2], v1, sep=":"),collapse=" + ")  
  f4 = formula(paste(paste(outcome, f1, sep=" ~ "),f3,sep=" + "))
  
  olsh = lm(f4, data)
  olsh$coef[is.na(olsh$coef)] = 0
  ateolsh1 = coef(olsh)[names(coef(olsh))==as.character(formula[2])]
  
  mat1 = model.matrix(formula(paste("~",paste(v1, collapse="+"),sep="")), data)
  mat1 = mat1[,colnames(mat1)!=c("(Intercept)")]
  ateolsh2 = colMeans(mat1) %*% olsh$coef[paste(as.character(formula)[2],colnames(mat1),sep=":")]
  ateolsh = ateolsh1 + ateolsh2
    
  # probit
  p1 = glm(formula, family = binomial(link = "probit"), data)
  data$p1 = pnorm(predict(p1, newdata=data))
  data = data[!is.na(data$p1),]
  # generate mills ratios
  data$sc1 = ifelse(data[,as.character(formula)[2]]==1,
                    -dnorm(predict(p1, newdata=data))/(pnorm(predict(p1, newdata=data))),
                    dnorm(predict(p1, newdata=data))/(1-pnorm(predict(p1, newdata=data))))
  data$sc11 = data[,as.character(formula)[2]]*data$sc1
  data$sc01 = (1-data[,as.character(formula)[2]])*data$sc1
    
  f5 = formula(paste(paste(outcome, f1, sep=" ~ ")," + sc11 + sc01",sep=""))
  # bvn regression
  ls = lm(f5, data) 
  # extract estimated rhos
  rs0 = coef(ls)[names(coef(ls))=="sc01"]
  rsd = coef(ls)[names(coef(ls))=="sc11"] - coef(ls)[names(coef(ls))=="sc01"]
  atebvn = coef(ls)[names(coef(ls))==as.character(formula)[2]]
  
  f6 = paste(paste(paste(outcome, f1, sep=" ~ "),f3,sep=" + ")," + sc11 + sc01",sep="")
  lsh = lm(f6, data) 
  # extract estimated rhos
  rs0h = coef(lsh)[names(coef(lsh))=="sc01"]
  rsdh = coef(lsh)[names(coef(lsh))=="sc11"]-coef(lsh)[names(coef(lsh))=="sc01"]
  lsh$coef[is.na(lsh$coef)] = 0
  atebvnh1 = lsh$coef[names(lsh$coef)==as.character(formula)[2]] 
  atebvnh2 = colMeans(mat1) %*% lsh$coef[paste(as.character(formula)[2],colnames(mat1),sep=":")]
    
  atebvnh = atebvnh1 + atebvnh2

  # perform grid search to min bias 
  grid = data.frame(h = seq(from = -5, to = 5, length.out = 1000)) 
  val1 = apply(grid, 1, function(x) abs((rs0 +(1-pnorm(x)) *rsd)*(dnorm(x)/(pnorm(x)*(1-pnorm(x))))))
  val2 = apply(grid, 1, function(x) abs((rs0h+(1-pnorm(x))*rsdh)*(dnorm(x)/(pnorm(x)*(1-pnorm(x))))))
  grid$val1 = val1 ; grid$val2 = val2
  # calculate pstar
  pstar = pnorm(grid[grid[,2]==min(grid[,2]),1])
  pstarh = pnorm(grid[grid[,3]==min(grid[,3]),1])
  
  # perform alpha search
  # for loop over values ranging from 0% to 100%
  # at each iteration loop records if alpha amount of...
  # treatment and not treatment groups in trimmed data
  alpha = cbind(seq(0,1,0.01),NA,NA)
  data$t1 = ifelse(data[as.character(formula)[2]]==1, 1, 0)
  data$nt1 = ifelse(data$t1==1, 0, 1)
  for(i in 1:dim(alpha)[1]){
    pcut = c(max(c(0.02,pstar-alpha[i,1])),min(c(0.98,pstar+alpha[i,1])))
    dattemp = data[data$p1>pcut[1] & data$p1<pcut[2],]
    pin = c(sum(dattemp$t1)/sum(data$t1), sum(dattemp$nt1)/sum(data$nt1))
    alpha[i,2] = ifelse(pin[1]>=theta & pin[2]>=theta,1,0)  
    
    pcut = c(max(c(0.02,pstarh-alpha[i,1])),min(c(0.98,pstarh+alpha[i,1])))
    dattemp = data[data$p1>pcut[1] & data$p1<pcut[2],]
    pin = c(sum(dattemp$t1)/sum(data$t1), sum(dattemp$nt1)/sum(data$nt1))
    alpha[i,3] = ifelse(pin[1]>=theta & pin[2]>=theta,1,0)
  }
  # drop alphas that dont satisfy requirments
  alpha1 = alpha[alpha[,2]==1,]
  alphah = alpha[alpha[,3]==1,]
  # what is the smallest alpha that does?
  alpha = min(alpha1[,1])
  alphah = min(alphah[,1])
  
  #cut points
  pcut = c(max(c(0.02,pstar-alpha)),min(c(0.98,pstar+alpha)))
  pcuth = c(max(c(0.02,pstarh-alphah)),min(c(0.98,pstarh+alphah)))
  
  datc = data.frame(data[data$p1>pcut[1] & data$p1<pcut[2],])
  datch = data.frame(data[data$p1>pcuth[1] & data$p1<pcuth[2],])
  
  # IPW ate trim above and below by 2%
  datipw = data[data$p1<=0.98 & data$p1>=0.02, ]
  datipw$d1 = (datipw[outcome]*datipw$t1)/datipw$p1
  datipw$d2 = datipw$t1/datipw$p1
  datipw$d3 = (datipw[outcome]*datipw$nt1)/(1-datipw$p1)
  datipw$d4 = (datipw$nt1)/(1-datipw$p1)
  ateipw = sum(datipw$d1)/sum(datipw$d2) - sum(datipw$d3)/sum(datipw$d4)
  
  # MB ate
  datc$d1 = (datc[outcome]*datc$t1)/datc$p1
  datc$d2 = datc$t1/datc$p1
  datc$d3 = (datc[outcome]*datc$nt1)/(1-datc$p1)
  datc$d4 = (datc$nt1)/(1-datc$p1)
  atemb = sum(datc$d1)/sum(datc$d2) - sum(datc$d3)/sum(datc$d4)
  
  datch$d1 = (datch[outcome]*datch$t1)/datch$p1
  datch$d2 = datch$t1/datch$p1
  datch$d3 = (datch[outcome]*datch$nt1)/(1-datch$p1)
  datch$d4 = (datch$nt1)/(1-datch$p1)
  atembh = sum(datch$d1)/sum(datch$d2) - sum(datch$d3)/sum(datch$d4)
  
  # BC ate
  Bate = -(rs0+(1-pstar)*rsd)*(dnorm(qnorm(pstar))/(pstar*(1-pstar)))
  atembbc = atemb-Bate
  
  Bateh = -(rs0h+(1-pstarh)*rsdh)*(dnorm(qnorm(pstarh))/(pstarh*(1-pstarh)))
  atembbch = atembh-Bateh
  
  p1a = ifelse(pnorm(predict(p1, newdata=data))==1,
               qnorm(0.99999), 
               predict(p1, newdata=data))  
  Bate1 = mean(-(rs0+(1-pnorm(p1a))*rsd) *
                  (dnorm(p1a)/(pnorm(p1a)*(1-pnorm(p1a)))))
  Bate2 = mean(-(rs0h+(1-pnorm(p1a))*rsdh) *
                  (dnorm(p1a)/(pnorm(p1a)*(1-pnorm(p1a)))))
  
  atebcipw = ateipw - Bate1
  atebcipwh = ateipw - Bate2

  cate = data.frame(cate = c(ateols,ateipw,atemb,atebvn,atembbc,atebcipw,
                     pstar,theta,pcut,dim(data)[1],dim(datc)[1]))

  hate = data.frame(hate = c(ateolsh,ateipw,atembh,atebvnh,atembbch,atebcipwh,
                           pstarh,theta,pcuth,dim(data)[1],dim(datch)[1]))
  
  rownames(cate) = rownames(hate) = c("OLS","IPW","MB","BVN","MB-BC","BC-IPW",
                                    "P*","Theta","P low","P high","N","N Trim")
  
  cate$cate = as.numeric(sprintf(paste("%.",paste(digits,"f",sep=""),sep=""),cate$cate))
  hate$hate = as.numeric(sprintf(paste("%.",paste(digits,"f",sep=""),sep=""),hate$hate))
  
  results = list(cate=cate, hate=hate, tprobit=p1)
  
  return(results)
}

bmate = function(x, ...) UseMethod("bmate")

bmate.default = function(formula, data, theta=0.25, digits=3, consate=TRUE){
  est = bmest(formula, data, theta, digits)
  est$call = match.call() 
  class(est) = "bmate"
  est$atetype = ifelse(consate==TRUE,"Constant","Heterogenous")
  est
}

print.bmate = function(x, ...){
  cat("Call:\n")
  print(x$call)

  
  if(x$atetype=="Constant"){
    cat("\n----------------------------------------------\n")  
    cat("Constant Average Treatment Effect Estimators")
    cat("\n----------------------------------------------\n")  
    cat("\nAverage Treatment Effects:\n")
    aa = x$cate[1:6,]
    names(aa) = c("OLS","IPW","MB","BVN","MB-BC","BC-IPW")
    print(aa)
    cat("\nTheta: ", x$cate[8,])  
    cat("\nBias Minimizing Propensity Score:", x$cate[7,])
    cat("\nLower and upper bounds for P: (",x$cate[9,],",",x$cate[10,],")")
    cat("\nInitial Sample Size:", x$cate[11,]) 
    cat("\nTrimmed Sample Size:", x$cate[12,])  
  } else {
    cat("\n----------------------------------------------\n")  
    cat("Heterogenous Average Treatment Effect Estimators")
    cat("\n----------------------------------------------\n")  
    cat("\nAverage Treatment Effects:\n")
    aa = x$hate[1:6,]
    names(aa) = c("OLS","IPW","MB","BVN","MB-BC","BC-IPW")
    print(aa)
    cat("\nTheta: ", x$hate[8,])  
    cat("\nBias Minimizing Propensity Score:", x$hate[7,])
    cat("\nLower and upper bounds for P: (",x$hate[9,],",",x$hate[10,],")")
    cat("\nInitial Sample Size:", x$hate[11,]) 
    cat("\nTrimmed Sample Size:", x$hate[12,])  
  }
}



  




