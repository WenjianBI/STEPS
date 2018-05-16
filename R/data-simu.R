#' @title Data Simulation
#' @description A function to simulate dataset of primary trait, secondary trait, genotype and one covariate
#' @param par.ls an R list of b0,b1,b3,g0,g1. More information can be seen in Details.
#' @param sec.type an R character to specify secodnary trait type: "binary" or "continuous"
#' @param sd1 error term standard deviation for primary trait
#' @param sd2 error term standard deviation for secondary trait
#' @param N sample size of dataset
#' @param maf minor allele frequency of SNPs to simulate genotype
#' @param cutoff cutoff to generate binary secondary phenotype
#' @param qntl quantile to choose y1 and y2, between 0 and 0.5. Default value is 0.1, that is, subjects with primary phenotype of top 10\% and bottom 10\% are in cohort
#' @return An R matrix with each row for one subject. Columns contain the following components:
#' 'Y' is for primary traits, 'Z'/'D' is for continuous/binary secondary traits, 'G' is for genotypes, 'E' is for covariates.
#' @examples
#' par.ls = list(b0=0,b1=rnorm(2),b3=rnorm(1),g0=0,g1=rnorm(2))
#' data.cont = data.simu(par.ls,"continuous")
#' data.bina = data.simu(par.ls,"binary")
#' @details
#' For continuous secondary traits, dataset is simulated by the following model
#' \deqn{Z = g0+g1[1]G+g1[2]X+e1}
#' \deqn{Y = b0+b1[1]G+b1[2]X+b3Z+e2}
#' For binary secondary traits, dataset is simulated by the following model
#' \deqn{Z = g0+g1[1]G+g1[2]X+e1}
#' \deqn{D = I(Z>cutoff)}
#' \deqn{Y = b0+b1[1]G+b1[2]X+b3Z+e2}
#' where 'Z'/'D' is continuous/binary secondary trait, 'Y' is primary trait, 'X' is covariate following standard normal distribution,
#' 'G' is genotype following HWE with MAF of 'maf', error term 'e1'/'e2' follows normal distribution with a mean of 0 and standard deviation of 'sd1'/'sd2',
#' only subjects with primary phenotype at top/bottom quantile of 'qntl' are retained as extreme phenotype sampling design.
#' @export
data.simu <- function(par.ls,
                      sec.type,
                      sd1=1,
                      sd2=1,
                      N=1000,
                      maf=0.3,
                      cutoff=0,
                      qntl=0.1)
{
  if(qntl>0.5) stop("Argument 'qntl' should be less than 0.5.")
  if(!all(is.element(c("b0","b1","b3","g0","g1"),names(par.ls)))) stop("Argument 'par.ls' should have elements of b0,b1,b3,g0,g1.")
  if(sec.type!="binary"&sec.type!="continuous") stop("Argument 'sec.type' should be either 'binary' or 'continuous'.")
  if(sec.type=="binary") data.out = data.simu.bina(par.ls,cutoff,sd1,N,maf,qntl)
  if(sec.type=="continuous") data.out = data.simu.cont(par.ls,sd1,sd2,N,maf,qntl)
  return(data.out)
}


#####################################################################
#
# -------- lower functions for continuous secondary traits ----------
#
#####################################################################


data.simu.cont <- function(par.ls,  # list of b0,b1,b3,g0,g1,
                           sd1=1,sd2=1,N=1000,maf=0.3,
                           qntl=0.1)       # quantile to choose y1 and y2
{
  # G is for genotype; E is for environment; e1(e2) is for error following normal distribution

  g0 = par.ls$g0; g1 = par.ls$g1;
  b0 = par.ls$b0; b1 = par.ls$b1; b3 = par.ls$b3;

  n.E = length(g1)-1  # number of covariates (envrionment)
  n = floor(N/(qntl*2));            # number of subjects
  G = rbinom(n,2,maf)
  E = matrix(rnorm(n*n.E),n,n.E)
  GE = cbind(G,E)
  err1 = rnorm(n,0,sd1)
  err2 = rnorm(n,0,sd2)

  # generate Y (primary trait) and Z (secondary trait)

  Z = g0+c(GE%*%g1)+err1
  Y = b0+c(GE%*%b1)+b3*Z+err2
  data = data.frame(Y=Y,Z=Z,G=G,E=E)
  rm(G,E,Z,Y)

  # y1 and y2 are cutoffs for extreme value sampling process

  q1 = qntl;   # quantile 1 for extremely low
  q2 = 1-qntl;   # quantile 2 for extremely high
  y1 = quantile(data$Y,q1)
  y2 = quantile(data$Y,q2)

  # extreme value sampling design

  ok = which(data$Y>y2|data$Y<y1)
  data.mat = as.matrix(data[ok,])
  rownames(data.mat) = NULL
  return(data.mat)
}
# example
# par.ls = list(b0=0,b1=rnorm(2),b3=rnorm(1),g0=0,g1=rnorm(2))
# data.mat = data.simu(par.ls)



#################################################################
#
# -------- lower functions for binary secondary traits ----------
#
#################################################################

data.simu.bina <- function(par.ls,cutoff,  # list of b0,b1,b3,g0,g1,
                           sd1=1,N=1000,maf=0.3,
                           qntl=0.1)       # quantile to choose y1 and y2
{
  # G is for genotype; E is for environment; e1(e2) is for error following normal distribution

  g0 = par.ls$g0; g1 = par.ls$g1;
  b0 = par.ls$b0; b1 = par.ls$b1; b3 = par.ls$b3;

  n.E = length(g1)-1  # number of covariates (envrionment)
  n = floor(N/(qntl*2));            # number of subjects
  G = rbinom(n,2,maf)
  E = matrix(rnorm(n*n.E),n,n.E)
  GE = cbind(G,E)
  err1 = rnorm(n,0,sd1)
  err2 = rnorm(n,0,1)

  # generate Y (primary trait) and Z (secondary trait)

  Z = g0+c(GE%*%g1)+err1
  D = as.numeric(Z>cutoff)
  Y = b0+c(GE%*%b1)+b3*D+err2
  data = data.frame(Y=Y,D=D,G=G,E=E)
  rm(G,E,Z,Y)

  # y1 and y2 are cutoffs for extreme value sampling process

  q1 = qntl;   # quantile 1 for extremely low
  q2 = 1-qntl;   # quantile 2 for extremely high
  y1 = quantile(data$Y,q1)
  y2 = quantile(data$Y,q2)

  # extreme value sampling design

  ok = which(data$Y>y2|data$Y<y1)
  data.mat = as.matrix(data[ok,])
  rownames(data.mat) = NULL
  return(data.mat)
}

