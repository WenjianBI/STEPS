
##################################################################################################################################
#
# -------- lower functions to compute log-likelihood and derivation over all parameters for continuous secondary traits ----------
#
##################################################################################################################################

ll.comp.cont = function(data.mat,   # numeric matrix of data with one row per subject. Column 1 is primary phneotype, column 2 is secondary phenotype, the other columns are genotype or covariates
                        y0=-Inf,    # lower cutoff less than y1.
                        y1,         # lower cutoff for primary phenotype in extreme sampling design
                        y2,         # upper cutoff for primary phenotype in extreme sampling design
                        y3=Inf,     # upper cufoff greater than y2.
                        para)       # a list of parameters containing (g0,g1,b0,b1,b3,sd1,sd2), note that g1 and b1 may be vector whose dimension is the same with the sum of genotypes and covariates
{
  # extract data
  Y = data.mat[,1]                 # primary phenotype
  Z = data.mat[,2]                 # secondary phenotype
  n.GE = ncol(data.mat)-2          # number of genotype and covariates
  GE = data.mat[,2+1:n.GE,drop=F]  # keep matrix form
  # extract parameter
  g0 = para$g0; g1 = para$g1;
  b0 = para$b0; b1 = para$b1; b3 = para$b3;
  sd1 = para$sd1; sd2 = para$sd2;
  if(sd1<0|sd2<0) return(list(ll=Inf,d.ll=NA))
  if(length(g1)!=n.GE|length(b1)!=n.GE) stop("Numbers of parameters and data are not matched. Please check the dimension.")
  # intermediate parameter
  t1 = b0+c(GE%*%b1);                     # numeric vector: theta1 = beta0+beta1*g+beta2*e+...
  t2 = g0+c(GE%*%g1);                     # numeric vector: theta2 = gamma0+gamma1*g+gamma2*e+...
  t3 = sqrt(b3^2*sd2^2+sd1^2)             # number: theta3 as new sd
  # an important equation: dnorm(Y-t1-b3*Z,sd=sd1)*dnorm(Z-t2)=dnorm(Y-t1-b3*t2,sd=t3)*dnorm(Z-mu(y),sd=sd) where mu(y)=(t2*sd1^2+(y-t1)*b3*sd2^2)/t3^2 and sd=sqrt(sd1^2*sd2^2/t3^2).
  t = t1+t2*b3                            # numeric vector: theta = theta1+theta2*beta3
  r1 = Y-t1-b3*Z;                         # numeric vector: rho1 = Y-theta1-beta3*Z
  r2 = Z-t2;                              # numeric vector: rho2 = Z-theta2
  r3 = y1-t;r3e=y0-t;
  r4 = y2-t;r4e=y3-t;
  ## Selection Probability: Pr(S=1|G,E)
  prob.s1=pnorm(r3,sd=t3)-pnorm(r3e,sd=t3)+pnorm(r4e,sd=t3)-pnorm(r4,sd=t3)
  # r5 is the derivative of log(prob.s1) over parameter theta3
  # r5 = (r3*dnorm(r3,sd=t3)-r4*dnorm(r4,sd=t3))/(pnorm(r3,sd=t3)+1-pnorm(r4,sd=t3))
  r5 = (r3*dnorm(r3,sd=t3)-r3e*dnorm(r3e,sd=t3)+r4e*dnorm(r4e,sd=t3)-r4*dnorm(r4,sd=t3))/prob.s1
  # r6 is the derivative of log(prob.s1) over parameter other than theta3
  # r6 = (dnorm(r3,sd=t3)-dnorm(r4,sd=t3))/(pnorm(r3,sd=t3)+1-pnorm(r4,sd=t3))
  r6 = (dnorm(r3,sd=t3)-dnorm(r3e,sd=t3)+dnorm(r4e,sd=t3)-dnorm(r4,sd=t3))/prob.s1
  # f is numeric vector for likelihood of all subjects
  # f = dnorm(r2,sd=sd2)*dnorm(r1,sd=sd1)/(pnorm(r3,sd=t3)+1-pnorm(r4,sd=t3))
  f = dnorm(r2,sd=sd2)*dnorm(r1,sd=sd1)/(pnorm(r3,sd=t3)-pnorm(r3e,sd=t3)+pnorm(r4e,sd=t3)-pnorm(r4,sd=t3))
  # derivative of f over parameter c(t1,t2,t3,b3,sd1,sd2)
  df.int = f*cbind(r1/sd1^2+r6,              # t1
                   r2/sd2^2+r6*b3,           # t2
                   r5/t3,                    # t3
                   r1/sd1^2*Z+r6*t2,         # b3
                   (r1^2-sd1^2)/sd1^3,       # sd1
                   (r2^2-sd2^2)/sd2^3)       # sd2
  a = cbind(1,GE)
  # derivative of f over all true parameters
  df = cbind(df.int[,2]*a,                             # g0,g1,g2,...
             df.int[,1]*a,                             # b0,b1,b2,...
             df.int[,3]*b3*sd2^2/t3+df.int[,4],        # b3
             df.int[,3]*sd1/t3+df.int[,5],             # sd1
             df.int[,3]*b3^2*sd2/t3+df.int[,6])        # sd2
  colnames(df) = c("g0",paste("g1",1:n.GE,sep="-"),
                   "b0",paste("b1",1:n.GE,sep="-"),
                   "b3","sd1","sd2")
  # mean of negative loglikelihood and corresponding derivative
  ll = mean(-1*log(f))
  d.ll = colMeans(-1*df/f)
  return(list(ll=ll,d.ll=d.ll))
}

##############################################################################################################################
#
# -------- lower functions to compute log-likelihood and derivation over all parameters for binary secondary traits ----------
#
##############################################################################################################################

ll.comp.bina = function(data.mat,   # numeric matrix of data with one row per subject. Column 1 is primary phneotype, column 2 is secondary phenotype, the other columns are genotype or covariates
                        y0=-Inf,    # lower cutoff less than y1.
                        y1,         # lower cutoff for primary phenotype in extreme sampling design
                        y2,         # upper cutoff for primary phenotype in extreme sampling design
                        y3=Inf,     # upper cufoff greater than y2.
                        para)       # a list of parameters containing (g0,g1,b0,b1,b3,sd1), note that g1 and b1 may be vector whose dimension is the same with the sum of genotypes and covariates
{
  # extract data
  Y = data.mat[,1]                 # primary phenotype
  D = data.mat[,2]                 # secondary phenotype
  n.GE = ncol(data.mat)-2          # number of genotype and covariates
  GE = data.mat[,2+1:n.GE,drop=F]  # keep matrix form
  # extract parameter
  g0 = para$g0; g1 = para$g1;
  b0 = para$b0; b1 = para$b1; b3 = para$b3;
  sd1 = para$sd1;
  if(sd1<0) return(list(ll=Inf,d.ll=NA))
  if(length(g1)!=n.GE|length(b1)!=n.GE) stop("Numbers of parameters and data are not matched. Please check the dimension.")
  # intermediate parameter
  t1 = b0+c(GE%*%b1);                     # numeric vector: theta1 = beta0+beta1*g+beta2*e+...
  t2 = g0+c(GE%*%g1);                     # numeric vector: theta2 = gamma0+gamma1*g+gamma2*e+...
  # an important equation: dnorm(Y-t1-b3*Z,sd=sd1)*dnorm(Z-t2)=dnorm(Y-t1-b3*t2,sd=t3)*dnorm(Z-mu(y),sd=sd) where mu(y)=(t2*sd1^2+(y-t1)*b3*sd2^2)/t3^2 and sd=sqrt(sd1^2*sd2^2/t3^2).
  r1 = pnorm(t2);
  r2 = dnorm(t2);
  r3 = Y-t1-b3*D;
  r3d = dnorm(r3,sd=sd1);
  r5 = cbind(y1-t1-b3,y1-t1,y2-t1-b3,y2-t1,y0-t1-b3,y0-t1,y3-t1-b3,y3-t1)
  colnames(r5)=c("y11","y12","y21","y22","y01","y02","y31","y32")
  r5d = dnorm(r5,sd=sd1)
  r5p = pnorm(r5,sd=sd1)
  ## A minor update to avoid Inf*0=NaN
  if(y0==-Inf) y0=-100
  if(y3==Inf) y3=100
  r5 = cbind(y1-t1-b3,y1-t1,y2-t1-b3,y2-t1,y0-t1-b3,y0-t1,y3-t1-b3,y3-t1)
  r6 = r5d*r5*(-1)/sd1;
  # ll is the mean of negative log-likelihood for all subjects
  f1 = r1*D+(1-r1)*(1-D)
  f2 = r3d
  f3 = r1*(r5p[,"y11"]-r5p[,"y01"]+r5p[,"y31"]-r5p[,"y21"])+(1-r1)*(r5p[,"y12"]-r5p[,"y02"]+r5p[,"y32"]-r5p[,"y22"])
  ll = mean(-1*log(f1*f2/f3))
  # derivative of f over parameter c(t1,t2,t3,b3,sd1,sd2)
  df1 = cbind(0,
              r2*(2*D-1),
              0,
              0)
  df2 = f2*cbind(r3/sd1^2,
                 0,
                 r3/sd1^2*D,
                 (r3^2-sd1^2)/sd1^3)
  df3 = cbind(-1*(r1*(r5d[,"y11"]-r5d[,"y01"]+r5d[,"y31"]-r5d[,"y21"])+(1-r1)*(r5d[,"y12"]-r5d[,"y02"]+r5d[,"y32"]-r5d[,"y22"])), # t1
              r2*(r5p[,"y11"]-r5p[,"y01"]+r5p[,"y31"]-r5p[,"y21"])-r2*(r5p[,"y12"]-r5p[,"y02"]+r5p[,"y32"]-r5p[,"y22"]),      # t2
              -1*(r1*(r5d[,"y11"]-r5d[,"y01"]+r5d[,"y31"]-r5d[,"y21"])),                          # b3
              r1*(r6[,"y11"]-r6[,"y01"]+r6[,"y31"]-r6[,"y21"])+(1-r1)*(r6[,"y12"]-r6[,"y02"]+r6[,"y32"]-r6[,"y22"]))          # sd1
  df.int = df3/f3-df1/f1-df2/f2
  a = cbind(1,GE)
  # derivative of f over all true parameters
  df = cbind(df.int[,2]*a,           # g0,g1,g2,...
             df.int[,1]*a,           # b0,b1,b2,...
             df.int[,3],             # b3
             df.int[,4])             # sd1
  colnames(df) = c("g0",paste("g1",1:n.GE,sep="-"),
                   "b0",paste("b1",1:n.GE,sep="-"),
                   "b3","sd1")
  # mean of negative loglikelihood and corresponding derivative
  d.ll = colMeans(df)
  return(list(ll=ll,d.ll=d.ll))
}
