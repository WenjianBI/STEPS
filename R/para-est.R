
###########################################################################################
#
# -------- lower functions to estimate parameter for continuous secondary traits ----------
#
###########################################################################################

para.est.cont = function(data.mat,       # numeric matrix of data with one row per subject. Column 1 is primary phneotype, column 2 is secondary phenotype, column 3 is genotype, the other columns are covariates
                         ini.para=NULL,  # a list of parameters containing (g0,g1,b0,b1,b3,sd1,sd2), note that g1 and b1 may be vector whose dimension is the same with the sum of genotypes and covariates
                         y0=-Inf,
                         y1=NULL,        # upper cutoff for primary phenotype in extreme sampling design
                         y2=NULL,        # lower cutoff for primary phenotype in extreme sampling design
                         y3=Inf,
                         method="BFGS")  # see function optim() for optional optimization method
{
  if(is.null(ini.para)) ini.para = ini.para.est.cont(data.mat);
  Y.c = sort(data.mat[,1])
  Y.c = Y.c[which(Y.c>quantile(Y.c,0.1)&Y.c<quantile(Y.c,0.9))]  # keep central 80% Y to compute cutoff for Y
  brk=which.max(diff(Y.c))
  if(is.null(y1)) y1 = Y.c[brk]          # lower cutoff
  if(is.null(y2)) y2 = Y.c[brk+1]        # upper cutoff
  if(y0==-Inf) y0=min(data.mat[,1])
  if(y3==Inf) y3=max(data.mat[,1])
  n.GE = ncol(data.mat)-2;
  par = c(ini.para$g0,unlist(ini.para$g1),ini.para$b0,unlist(ini.para$b1),ini.para$b3,ini.para$sd1,ini.para$sd2);
  fn = function(par.v){
    para = list(g0 = par.v[1], g1 = par.v[1+1:n.GE],
                b0 = par.v[2+n.GE], b1 = par.v[2+n.GE+1:n.GE], b3 = par.v[3+2*n.GE],
                sd1 = par.v[4+2*n.GE], sd2 = par.v[5+2*n.GE])
    res.ll = ll.comp.cont(data.mat,y0,y1,y2,y3,para)
    return(res.ll$ll)
  }

  gr = function(par.v){
    para = list(g0 = par.v[1], g1 = par.v[1+1:n.GE],
                b0 = par.v[2+n.GE], b1 = par.v[2+n.GE+1:n.GE], b3 = par.v[3+2*n.GE],
                sd1 = par.v[4+2*n.GE], sd2 = par.v[5+2*n.GE])
    res.ll = ll.comp.cont(data.mat,y0,y1,y2,y3,para)
    return(res.ll$d.ll)
  }

  res.opt=optim(par,fn,gr,method = method)
  names(res.opt$par) = c("g0",paste("g1",1:n.GE,sep="-"),
                         "b0",paste("b1",1:n.GE,sep="-"),
                         "b3","sd1","sd2")
  res.lst = list(g0=res.opt$par[1],g1=res.opt$par[1+1:n.GE],
                 b0=res.opt$par[2+n.GE],b1=res.opt$par[2+n.GE+1:n.GE],b3=res.opt$par[3+2*n.GE],
                 sd1=res.opt$par[4+2*n.GE],sd2=res.opt$par[5+2*n.GE])
  return(list(res.opt=res.opt,res.lst=res.lst))
}

#######################################################################################
#
# -------- lower functions to estimate parameter for binary secondary traits ----------
#
#######################################################################################

para.est.bina = function(data.mat,       # numeric matrix of data with one row per subject. Column 1 is primary phneotype, column 2 is secondary phenotype, column 3 is genotype, the other columns are covariates
                         ini.para=NULL,  # a list of parameters containing (g0,g1,b0,b1,b3,sd1), note that g1 and b1 may be vector whose dimension is the same with the sum of genotypes and covariates
                         y0=-Inf,
                         y1=NULL,        # upper cutoff for primary phenotype in extreme sampling design
                         y2=NULL,        # lower cutoff for primary phenotype in extreme sampling design
                         y3=Inf,
                         method="BFGS")  # see function optim() for optional optimization method
{
  if(is.null(ini.para)) ini.para = ini.para.est.bina(data.mat);
  Y.c = sort(data.mat[,1])
  Y.c = Y.c[which(Y.c>quantile(Y.c,0.1)&Y.c<quantile(Y.c,0.9))]  # keep central 80% Y to compute cutoff for Y
  brk=which.max(diff(Y.c))
  if(is.null(y1)) y1 = Y.c[brk]
  if(is.null(y2)) y2 = Y.c[brk+1]
  n.GE = ncol(data.mat)-2;
  par = c(ini.para$g0,unlist(ini.para$g1),ini.para$b0,unlist(ini.para$b1),ini.para$b3,ini.para$sd1);
  fn = function(par.v){
    para = list(g0 = par.v[1], g1 = par.v[1+1:n.GE],
                b0 = par.v[2+n.GE], b1 = par.v[2+n.GE+1:n.GE], b3 = par.v[3+2*n.GE],
                sd1 = par.v[4+2*n.GE])
    res.ll = ll.comp.bina(data.mat,y0,y1,y2,y3,para)
    return(res.ll$ll)
  }

  gr = function(par.v){
    para = list(g0 = par.v[1], g1 = par.v[1+1:n.GE],
                b0 = par.v[2+n.GE], b1 = par.v[2+n.GE+1:n.GE], b3 = par.v[3+2*n.GE],
                sd1 = par.v[4+2*n.GE])
    res.ll = ll.comp.bina(data.mat,y0,y1,y2,y3,para)
    return(res.ll$d.ll)
  }

  res.opt=optim(par,fn,gr,method = method,control = list(maxit=1000,reltol=1e-32))
  names(res.opt$par) = c("g0",paste("g1",1:n.GE,sep="-"),
                         "b0",paste("b1",1:n.GE,sep="-"),
                         "b3","sd1")
  res.lst = list(g0=res.opt$par[1],g1=res.opt$par[1+1:n.GE],
                 b0=res.opt$par[2+n.GE],b1=res.opt$par[2+n.GE+1:n.GE],b3=res.opt$par[3+2*n.GE],
                 sd1=res.opt$par[4+2*n.GE])
  return(list(res.opt=res.opt,res.lst=res.lst))
}


# lower function. Directly give a simple estimate of parameter as the initial parameter for the further iterative algorithm
ini.para.est.cont = function(data.mat)
{
  Y = data.mat[,1]  # primary phenotype
  Z = data.mat[,2]  # secondary phenotype
  n.GE = ncol(data.mat)-2  # number of genotype and covariates
  GE = data.mat[,2+1:n.GE,drop=F]  # keep matrix form
  res.Y = lm(Y~GE+Z)
  b0 = res.Y$coefficients[1]
  b1 = res.Y$coefficients[1+1:n.GE]
  b3 = res.Y$coefficients[2+n.GE]
  sd1 = summary(res.Y)$sigma
  res.Z = lm(Z~GE)
  g0 = res.Z$coefficients[1]
  g1 = res.Z$coefficients[1+1:n.GE]
  sd2 = summary(res.Z)$sigma
  return(list(g0=g0,g1=g1,b0=b0,b1=b1,b3=b3,sd1=sd1,sd2=sd2))
}


# lower function. Directly give a simple estimate of parameter as the initial parameter for the further iterative algorithm
ini.para.est.bina = function(data.mat)
{
  #   Y = data.mat[,1]  # primary phenotype
  #   D = data.mat[,2]  # secondary phenotype
  #   n.GE = ncol(data.mat)-2  # number of genotype and covariates
  #   GE = data.mat[,2+1:n.GE,drop=F]  # keep matrix form
  #   res.Y = lm(Y~GE+D)
  #   b0 = res.Y$coefficients[1]
  #   b1 = res.Y$coefficients[1+1:n.GE]
  #   b3 = res.Y$coefficients[2+n.GE]
  #   sd1 = summary(res.Y)$sigma
  #   res.Z = glm(D~GE,family = binomial)
  #   g0 = res.Z$coefficients[1]
  #   g1 = res.Z$coefficients[1+1:n.GE]
  #   return(list(g0=g0,g1=g1,b0=b0,b1=b1,b3=b3,sd1=sd1))
  n.GE = ncol(data.mat)-2  # number of genotype and covariates
  return(list(g0=0,g1=rep(0,n.GE),b0=0,b1=rep(0,n.GE),b3=0,sd1=1))
}






