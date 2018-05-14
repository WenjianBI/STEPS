
##############################################################################################
#
# -------- lower functions to compute fisher matrix for continuous secondary traits ----------
#
##############################################################################################

fisher.comp.cont <- function(data.mat,   # data in matrix where each row is for one observation, first three columns are for primary phneotype, secondary phenotype, and genotype, the other columns are for covariates
                             y0=-Inf,
                             y1=NULL,         # upper cutoff for primary phenotype in extreme sampling design
                             y2=NULL,         # lower cutoff for primary phenotype in extreme sampling design
                             y3=Inf,
                             para)       # a list of parameters containing (g0,g1,b0,b1,b3,sd1,sd2), note that g1 and b1 may be vector whose dimension is the same with the sum of genotypes and covariates
{
  Y.c = sort(data.mat[,1])
  Y.c = Y.c[which(Y.c>quantile(Y.c,0.1)&Y.c<quantile(Y.c,0.9))]  # only keep central 80% Y
  brk=which.max(diff(Y.c))
  y1.est = Y.c[brk]
  y2.est = Y.c[brk+1]
  if(is.null(y1)) y1=y1.est
  if(is.null(y2)) y2=y2.est
  n.GE = ncol(data.mat)-2
  GE = data.mat[,2+1:n.GE,drop=F]  # genotype and covariate
  # parameter
  g0 = para$g0; g1 = para$g1;
  b0 = para$b0; b1 = para$b1; b3 = para$b3;
  sd1 = para$sd1; sd2 = para$sd2;
  # theta1, theta2, and theta3
  t1 = b0+c(GE%*%b1);                     # theta1 = beta0+beta1*g+beta2*e
  t2 = g0+c(GE%*%g1);                     # theta2 = gamma0+gamma1*g+gamma2*e
  t3 = sqrt(b3^2*sd2^2+sd1^2)                 # theta3 : number
  # pdf and cdf for some points
  t = t1+t2*b3
  F1 = pnorm(y1-t,sd=t3);F1e=pnorm(y0-t,sd=t3)
  F2 = pnorm(y2-t,sd=t3);F2e=pnorm(y3-t,sd=t3)
  f1 = dnorm(y1-t,sd=t3);f1e=dnorm(y0-t,sd=t3)
  f2 = dnorm(y2-t,sd=t3);f2e=dnorm(y3-t,sd=t3)
  ## A minor update to avoid Inf*0=NaN
  if(y0==-Inf) y0=-100
  if(y3==Inf) y3=100
  # prop=1/(F1+1-F2)
  prop = 1/(F1-F1e+F2e-F2);   # proportion of 1/(F1+1-F2), a vector
  # mu=mu1*y+mu2,
  mu1 = (b3*sd2^2)/t3^2
  mu2 = (t2*sd1^2-t1*b3*sd2^2)/t3^2
  sd = sd1*sd2/abs(t3)
  # matrix of 15 by N to show integrate of y*p(..)
  N = nrow(data.mat)
  int.mat = matrix(NA,15,N,dimnames = list(c("1","y","y^2","y^3","y^4",   # 1-5
                                             "z","z*y","z*y^2","z*y^3",   # 6-9
                                             "z^2","z^2*y","z^2*y^2",     # 10-12
                                             "z^3","z^3*y","z^4"),NULL))  # 13-15
  # int.mat[1,] = 1-F2+F1;                                                  #1        1
  int.mat[1,] = F1-F1e+F2e-F2
  # int.mat[2,] = t*int.mat[1,]-t3^2*(f1-f2)                                #y        2
  int.mat[2,] = t*int.mat[1,]-t3^2*(f1-f1e+f2e-f2)
  # int.mat[3,] = t*int.mat[2,]+t3^2*int.mat[1,]-t3^2*(f1*y1-f2*y2)         #y^2      3
  int.mat[3,] = t*int.mat[2,]+t3^2*int.mat[1,]-t3^2*(f1*y1-f1e*y0+f2e*y3-f2*y2)
  # int.mat[4,] = t*int.mat[3,]+2*t3^2*int.mat[2,]-t3^2*(y1^2*f1-y2^2*f2)   #y^3      4
  int.mat[4,] = t*int.mat[3,]+2*t3^2*int.mat[2,]-t3^2*(f1*y1^2-f1e*y0^2+f2e*y3^2-f2*y2^2)
  # int.mat[5,] = t*int.mat[4,]+3*t3^2*int.mat[3,]-t3^2*(y1^3*f1-y2^3*f2)   #y^4      5
  int.mat[5,] = t*int.mat[4,]+3*t3^2*int.mat[3,]-t3^2*(f1*y1^3-f1e*y0^3+f2e*y3^3-f2*y2^3)

  int.mat[6,] = mu1*int.mat[2,]+mu2*int.mat[1,]                           #z = mu1*y+mu2
  int.mat[7,] = mu1*int.mat[3,]+mu2*int.mat[2,]                           #z*y = (mu1*y+mu2)*y=mu1*y^2+mu2*y
  int.mat[8,] = mu1*int.mat[4,]+mu2*int.mat[3,]                           #z*y^2 = mu1*y^3+mu2*y^2
  int.mat[9,] = mu1*int.mat[5,]+mu2*int.mat[4,]                           #z*y^3 = mu1*y^4+mu2*y^3
  int.mat[10,] = mu1^2*int.mat[3,]+(2*mu1*mu2)*int.mat[2,]+(mu2^2+sd^2)*int.mat[1,]  #z^2 = (mu1*y+mu2)^2+sd^2=mu1^2*y^2+(2*mu1*mu2)*y+(mu2^2+sd^2)
  int.mat[11,] = mu1^2*int.mat[4,]+(2*mu1*mu2)*int.mat[3,]+(mu2^2+sd^2)*int.mat[2,]  #z^2*y = mu1^2*y^3+(2*mu1*mu2)*y^2+(mu2^2+sd^2)*y
  int.mat[12,] = mu1^2*int.mat[5,]+(2*mu1*mu2)*int.mat[4,]+(mu2^2+sd^2)*int.mat[3,]  #z^2*y^2 = mu1^2*y^4+(2*mu1*mu2)*y^3+(mu2^2+sd^2)*y^2
  int.mat[13,] = mu1^3*int.mat[4,]+(3*mu1^2*mu2)*int.mat[3,]+(3*mu1*mu2^2+3*mu1*sd^2)*int.mat[2,]+(mu2^3+3*mu2*sd^2)*int.mat[1,]  #z^3 = (mu1*y+mu2)^3+3*(mu1*y+mu2)*sd^2=mu1^3*y^3+(3*mu1^2*mu2)*y^2+(3*mu1*mu2^2+3*mu1*sd^2)*y+(mu2^3+3*mu2*sd^2)
  int.mat[14,] = mu1^3*int.mat[5,]+(3*mu1^2*mu2)*int.mat[4,]+(3*mu1*mu2^2+3*mu1*sd^2)*int.mat[3,]+(mu2^3+3*mu2*sd^2)*int.mat[2,]  #z^3*y = mu1^3*y^4+(3*mu1^2*mu2)*y^3+(3*mu1*mu2^2+3*mu1*sd^2)*y^2+(mu2^3+3*mu2*sd^2)*y
  int.mat[15,] = mu1^4*int.mat[5,]+(4*mu2*mu1^3)*int.mat[4,]+(6*mu2^2*mu1^2+6*sd^2*mu1^2)*int.mat[3,]+(4*mu2^3*mu1+12*sd^2*mu1*mu2)*int.mat[2,]+(mu2^4+6*sd^2*mu2^2+3*sd^4)*int.mat[1,]  #z^4 = (mu1*y+mu2)^4+6*(mu1*y+mu2)^2*sd^2+3*sd^4=mu1^4*y^4+(4*mu2*mu1^3)*y^3+(6*mu2^2*mu1^2+6*sd^2*mu1^2)*y^2+(4*mu2^3*mu1+12*sd^2*mu1*mu2)*y+(mu2^4+12*sd^2*mu1*mu2+sd^4)
  # matrix of 6 by N
  r3 = y1-t;r3e=y0-t
  r4 = y2-t;r4e=y3-t
  # r5 = (r3*f1-r4*f2)/(F1+1-F2)
  r5 = (r3*f1-r3e*f1e+r4e*f2e-r4*f2)/(F1-F1e+F2e-F2)
  # r6 = (f1-f2)/(F1+1-F2)
  r6 = (f1-f1e+f2e-f2)/(F1-F1e+F2e-F2)
  rho.mat = matrix(NA,15,N,dimnames = list(c("1","r1","r1^2","r1^3","r1^4",     # 1-5
                                             "z","z*r1","z*r1^2","z*r1^3",      # 6-9
                                             "z^2","z^2*r1","z^2*r1^2",         # 10-12
                                             "z^3","z^3*r1","z^4"),NULL))       # 13-15
  rho.mat[1,] = int.mat["1",]
  rho.mat[2,] = (-1*t1)*int.mat["1",]+int.mat["y",]+(-1*b3)*int.mat["z",];
  rho.mat[3,] = (t1^2)*int.mat["1",]+(-2*t1)*int.mat["y",]+int.mat["y^2",]+
    (2*t1*b3)*int.mat["z",]+(-2*b3)*int.mat["z*y",]+
    (b3^2)*int.mat["z^2",];
  rho.mat[4,] = (-1*t1^3)*int.mat["1",]+(3*t1^2)*int.mat["y",]+(-3*t1)*int.mat["y^2",]+int.mat["y^3",]+
    (-3*t1^2*b3)*int.mat["z",]+(6*t1*b3)*int.mat["z*y",]+(-3*b3)*int.mat["z*y^2",]+
    (-3*t1*b3^2)*int.mat["z^2",]+(3*b3^2)*int.mat["z^2*y",]+
    (-1*b3^3)*int.mat["z^3",];
  rho.mat[5,] = (t1^4)*int.mat["1",]+(-4*t1^3)*int.mat["y",]+(6*t1^2)*int.mat["y^2",]+(-4*t1)*int.mat["y^3",]+int.mat["y^4",]+
    (4*t1^3*b3)*int.mat["z",]+(-12*t1^2*b3)*int.mat["z*y",]+(12*t1*b3)*int.mat["z*y^2",]+(-4*b3)*int.mat["z*y^3",]+
    (6*t1^2*b3^2)*int.mat["z^2",]+(-12*b3^2*t1)*int.mat["z^2*y",]+(6*b3^2)*int.mat["z^2*y^2",]+
    (4*t1*b3^3)*int.mat["z^3",]+(-4*b3^3)*int.mat["z^3*y",]+
    (b3^4)*int.mat["z^4",];
  rho.mat[6,] = int.mat["z",];
  rho.mat[7,] = -t1*int.mat["z",]+int.mat["z*y",]-b3*int.mat["z^2",];
  rho.mat[8,] = t1^2*int.mat["z",]-2*t1*int.mat["z*y",]+int.mat["z*y^2",]+
    2*t1*b3*int.mat["z^2",]-2*b3*int.mat["z^2*y",]+
    b3^2*int.mat["z^3",]
  rho.mat[9,] = (-1)*t1^3*int.mat["z",]+3*t1^2*int.mat["z*y",]-3*t1*int.mat["z*y^2",]+int.mat["z*y^3",]+
    (-3)*t1^2*b3*int.mat["z^2",]+6*t1*b3*int.mat["z^2*y",]-3*b3*int.mat["z^2*y^2",]+
    (-3)*t1*b3^2*int.mat["z^3",]+3*b3^2*int.mat["z^3*y",]+
    (-1)*b3^3*int.mat["z^4",];
  rho.mat[10,] = int.mat["z^2",]
  rho.mat[11,] = -t1*int.mat["z^2",]+int.mat["z^2*y",]-b3*int.mat["z^3",];
  rho.mat[12,] = t1^2*int.mat["z^2",]-2*t1*int.mat["z^2*y",]+int.mat["z^2*y^2",]+
    2*t1*b3*int.mat["z^3",]-2*b3*int.mat["z^3*y",]+
    b3^2*int.mat["z^4",]
  rho.mat[13,] = int.mat["z^3",]
  rho.mat[14,] = -t1*int.mat["z^3",]+int.mat["z^3*y",]-b3*int.mat["z^4",];
  rho.mat[15,] = int.mat["z^4",]

  a1 = 1/sd1^2; a2 = r6;
  bb1 = 1/sd2^2; bb2 = -1*t2/sd2^2+r6*b3;
  c1 = r5/t3
  d1 = 1/sd1^2; d2 = r6*t2
  e1 = 1/sd1^3; e2 = -1/sd1
  f1 = 1/sd2^3; f2 = -2*t2/sd2^3; f3 = (t2^2-sd2^2)/sd2^3;
  D = matrix(0,6,N,dimnames = list(c("t1","t2","t3","b3","sd1","sd2"),NULL))
  D[1,] = a1*rho.mat["r1",]+a2*rho.mat["1",];
  D[2,] = bb1*rho.mat["z",]+bb2*rho.mat["1",];
  D[3,] = c1*rho.mat["1",];
  D[4,] = d1*rho.mat["z*r1",]+d2*rho.mat["1",]
  D[5,] = e1*rho.mat["r1^2",]+e2*rho.mat["1",]
  D[6,] = f1*rho.mat["z^2",]+f2*rho.mat["z",]+f3*rho.mat["1",];
  DD = matrix(0,21,N,dimnames = list(c("t1*t1","t1*t2","t1*t3","t1*b3","t1*sd1","t1*sd2",   # 1-6
                                       "t2*t2","t2*t3","t2*b3","t2*sd1","t2*sd2",           # 7-11
                                       "t3*t3","t3*b3","t3*sd1","t3*sd2",                   # 12-15
                                       "b3*b3","b3*sd1","b3*sd2",                           # 16-18
                                       "sd1*sd1","sd1*sd2",                                 # 19,20
                                       "sd2*sd2"),NULL))                                    # 21

  DD[1,] = a1*(a1*rho.mat["r1^2",]+a2*rho.mat["r1",])+a2*(a1*rho.mat["r1",]+a2*rho.mat["1",]);
  DD[2,] = bb1*(a1*rho.mat["z*r1",]+a2*rho.mat["z",])+bb2*(a1*rho.mat["r1",]+a2*rho.mat["1",]);
  DD[3,] = c1*D[1,]
  DD[4,] = d1*(a1*rho.mat["z*r1^2",]+a2*rho.mat["z*r1",])+d2*(a1*rho.mat["r1",]+a2*rho.mat["1",]);
  DD[5,] = e1*a1*rho.mat["r1^3",]+e1*a2*rho.mat["r1^2",]+e2*a1*rho.mat["r1",]+e2*a2*rho.mat["1",];
  DD[6,] = f1*a1*rho.mat["z^2*r1",]+f1*a2*rho.mat["z^2",]+f2*a1*rho.mat["z*r1",]+f2*a2*rho.mat["z",]+f3*a1*rho.mat["r1",]+f3*a2*rho.mat["1",];
  DD[7,] = bb1*(bb1*rho.mat["z^2",]+bb2*rho.mat["z",])+bb2*(bb1*rho.mat["z",]+bb2*rho.mat["1",]);
  DD[8,] = c1*D[2,];
  DD[9,] = d1*(bb1*rho.mat["z^2*r1",]+bb2*rho.mat["z*r1",])+d2*(bb1*rho.mat["z",]+bb2*rho.mat["1",]);
  DD[10,] = e1*(bb1*rho.mat["z*r1^2",]+bb2*rho.mat["r1^2",])+e2*(bb1*rho.mat["z",]+bb2*rho.mat["1",]);
  DD[11,] = f1*(bb1*rho.mat["z^3",]+bb2*rho.mat["z^2",])+f2*(bb1*rho.mat["z^2",]+bb2*rho.mat["z",])+f3*(bb1*rho.mat["z",]+bb2*rho.mat["1",]);
  DD[12,] = c1*D[3,];
  DD[13,] = c1*D[4,];
  DD[14,] = c1*D[5,];
  DD[15,] = c1*D[6,];
  DD[16,] = d1*(d1*rho.mat["z^2*r1^2",]+d2*rho.mat["z*r1",])+d2*(d1*rho.mat["z*r1",]+d2*rho.mat["1",])
  DD[17,] = e1*(d1*rho.mat["z*r1^3",]+d2*rho.mat["r1^2",])+e2*(d1*rho.mat["z*r1",]+d2*rho.mat["1",])
  DD[18,] = f1*(d1*rho.mat["z^3*r1",]+d2*rho.mat["z^2",])+f2*(d1*rho.mat["z^2*r1",]+d2*rho.mat["z",])+f3*(d1*rho.mat["z*r1",]+d2*rho.mat["1",])
  DD[19,] = e1*(e1*rho.mat["r1^4",]+e2*rho.mat["r1^2",])+e2*(e1*rho.mat["r1^2",]+e2*rho.mat["1",])
  DD[20,] = f1*(e1*rho.mat["z^2*r1^2",]+e2*rho.mat["z^2",])+f2*(e1*rho.mat["z*r1^2",]+e2*rho.mat["z",])+f3*(e1*rho.mat["r1^2",]+e2*rho.mat["1",])
  DD[21,] = f1*(f1*rho.mat["z^4",]+f2*rho.mat["z^3",]+f3*rho.mat["z^2",])+f2*(f1*rho.mat["z^3",]+f2*rho.mat["z^2",]+f3*rho.mat["z",])+f3*(f1*rho.mat["z^2",]+f2*rho.mat["z",]+f3*rho.mat["1",])

  f = t(cbind(1,GE))
  nf = nrow(f)   # (number of covariates) + 2
  Df = rbind(matrix(D[2,],nf,N,byrow=T)*f,   # 9*N
             matrix(D[1,],nf,N,byrow=T)*f,
             D[3,]*b3*sd2^2/t3+D[4,],        # b3
             D[3,]*sd1/t3+D[5,],             # sd1
             D[3,]*b3^2*sd2/t3+D[6,])        # sd2
  Df = Df*matrix(prop,nrow(Df),ncol(Df),byrow=T)
  fisher.1a = Df %*% t(Df)
  rsum = rowSums(Df)
  fisher.1b = rsum %*% t(rsum)
  fisher.1 = fisher.1b - fisher.1a;
  fisher.2 = matrix(0,nrow(fisher.1),ncol(fisher.1))
  for(i in 1:N){
    ge = GE[i,];
    a = c(1,ge)
    A = rbind(cbind(0,a,0,0,0,0),
              cbind(a,0,0,0,0,0),
              c(rep(0,2),c(b3*sd2^2/t3,1),rep(0,2)),
              c(rep(0,2),sd1/t3,0,1,0),
              c(rep(0,2),b3^2*sd2/t3,rep(0,2),1))
    B = matrix(NA,6,6)
    B[lower.tri(B,diag=T)]=DD[,i];
    B[upper.tri(B)] = t(B)[upper.tri(B)]
    fisher.2 = fisher.2 + prop[i]*A%*%B%*%t(A)
  }
  fisher = fisher.1+fisher.2
  colnames(fisher) = rownames(fisher) = c("g0",paste("g1",1:n.GE,sep="-"),
                                          "b0",paste("b1",1:n.GE,sep="-"),
                                          "b3","sd1","sd2")
  return(fisher)
}

##########################################################################################
#
# -------- lower functions to compute fisher matrix for binary secondary traits ----------
#
##########################################################################################

fisher.comp.bina = function(data.mat,   # data in matrix where each row is for one observation, first three columns are for primary phneotype, secondary phenotype, and genotype, the other columns are for covariates
                            y0=-Inf,
                            y1=NULL,         # upper cutoff for primary phenotype in extreme sampling design
                            y2=NULL,         # lower cutoff for primary phenotype in extreme sampling design
                            y3=Inf,
                            para)       # a list of parameters containing (g0,g1,b0,b1,b3,sd1,sd2), note that g1 and b1 may be vector whose dimension is the same with the sum of genotypes and covariates
{
  Y.c = sort(data.mat[,1])
  Y.c = Y.c[which(Y.c>quantile(Y.c,0.1)&Y.c<quantile(Y.c,0.9))]  # only keep central 80% Y
  brk=which.max(diff(Y.c))
  y1.est = Y.c[brk]
  y2.est = Y.c[brk+1]
  if(is.null(y1)) y1=y1.est
  if(is.null(y2)) y2=y2.est
  n.GE = ncol(data.mat)-2
  GE = data.mat[,2+1:n.GE,drop=F]  # genotype and covariate
  # parameter
  g0 = para$g0; g1 = para$g1;
  b0 = para$b0; b1 = para$b1; b3 = para$b3;
  sd1 = para$sd1;
  # intermediate parameter
  t1 = b0+c(GE%*%b1);                     # numeric vector: theta1 = beta0+beta1*g+beta2*e+...
  t2 = g0+c(GE%*%g1);                     # numeric vector: theta2 = gamma0+gamma1*g+gamma2*e+...
  # an important equation: dnorm(Y-t1-b3*Z,sd=sd1)*dnorm(Z-t2)=dnorm(Y-t1-b3*t2,sd=t3)*dnorm(Z-mu(y),sd=sd) where mu(y)=(t2*sd1^2+(y-t1)*b3*sd2^2)/t3^2 and sd=sqrt(sd1^2*sd2^2/t3^2).
  r1 = pnorm(t2);
  r1s = pnorm(t2,lower.tail = F)  # 1-r1
  r2 = dnorm(t2);
  # r5 = cbind(y1-t1-b3,y1-t1,y2-t1-b3,y2-t1)
  r5 = cbind(y1-t1-b3,y1-t1,y2-t1-b3,y2-t1,y0-t1-b3,y0-t1,y3-t1-b3,y3-t1)
  colnames(r5)=c("y11","y12","y21","y22","y01","y02","y31","y32")
  r5d = dnorm(r5,sd=sd1)
  r5p = pnorm(r5,sd=sd1)
  f3 = r1*(r5p[,"y11"]-r5p[,"y01"]+r5p[,"y31"]-r5p[,"y21"])+(1-r1)*(r5p[,"y12"]-r5p[,"y02"]+r5p[,"y32"]-r5p[,"y22"])  # numeric vector: P(s=1)
  prop=1/f3

  ## A minor update to avoid Inf*0=NaN
  if(y0==-Inf) y0=-100
  if(y3==Inf) y3=100
  r5 = cbind(y1-t1-b3,y1-t1,y2-t1-b3,y2-t1,y0-t1-b3,y0-t1,y3-t1-b3,y3-t1)
  colnames(r5)=c("y11","y12","y21","y22","y01","y02","y31","y32")
  r6 = r5d*r5*(-1)/sd1;
  # ll is the mean of negative log-likelihood for all subjects
  N = nrow(data.mat)
  int.mat.1 = matrix(NA,N,5)  # d=1: N*5 matrix to integrate over y within (-Inf,y1)U(y2,Inf) for (y-t1-b3*d)^k where k=0,1,2,3,4
  int.mat.1[,1]=r5p[,"y31"]-r5p[,"y21"]+r5p[,"y11"]-r5p[,"y01"];
  int.mat.1[,2]=sd1^2*(r5d[,"y21"]-r5d[,"y31"]+r5d[,"y01"]-r5d[,"y11"]);
  int.mat.1[,3]=sd1^2*(int.mat.1[,1]+r5[,"y21"]*r5d[,"y21"]-r5[,"y31"]*r5d[,"y31"]+r5[,"y01"]*r5d[,"y01"]-r5[,"y11"]*r5d[,"y11"])
  int.mat.1[,4]=sd1^2*(2*int.mat.1[,2]+r5[,"y21"]^2*r5d[,"y21"]-r5[,"y31"]^2*r5d[,"y31"]+r5[,"y01"]^2*r5d[,"y01"]-r5[,"y11"]^2*r5d[,"y11"])
  int.mat.1[,5]=sd1^2*(3*int.mat.1[,3]+r5[,"y21"]^3*r5d[,"y21"]-r5[,"y31"]^3*r5d[,"y31"]+r5[,"y01"]^3*r5d[,"y01"]-r5[,"y11"]^3*r5d[,"y11"])

  int.mat.0 = matrix(NA,N,5)  # d=0: N*5 matrix to integrate over y within (-Inf,y1)U(y2,Inf) for (y-t1-b3*d)^k where k=0,1,2,3,4
  int.mat.0[,1]=r5p[,"y32"]-r5p[,"y22"]+r5p[,"y12"]-r5p[,"y02"]
  int.mat.0[,2]=sd1^2*(r5d[,"y22"]-r5d[,"y32"]+r5d[,"y02"]-r5d[,"y12"])
  int.mat.0[,3]=sd1^2*(int.mat.0[,1]+r5[,"y22"]*r5d[,"y22"]-r5[,"y32"]*r5d[,"y32"]+r5[,"y02"]*r5d[,"y02"]-r5[,"y12"]*r5d[,"y12"])
  int.mat.0[,4]=sd1^2*(2*int.mat.0[,2]+r5[,"y22"]^2*r5d[,"y22"]-r5[,"y32"]^2*r5d[,"y32"]+r5[,"y02"]^2*r5d[,"y02"]-r5[,"y12"]^2*r5d[,"y12"])
  int.mat.0[,5]=sd1^2*(3*int.mat.0[,3]+r5[,"y22"]^3*r5d[,"y22"]-r5[,"y32"]^3*r5d[,"y32"]+r5[,"y02"]^3*r5d[,"y02"]-r5[,"y12"]^3*r5d[,"y12"])

  int.mat = int.mat.1*r1+int.mat.0*(1-r1)  # N*5 matrix to integrate over d and y for (y-t1-b3*d)^k where k=0,1,2,3,4

  # N*4  integrate of derivative of log(P(y,d)) over y and d within (-Inf,y1)U(y2,Inf)
  A = cbind(int.mat[,2]/sd1^2,                                    # t1
            int.mat.1[,1]*r2+int.mat.0[,1]*(-1)*r2,               # t2
            int.mat.1[,2]/sd1^2*r1,                               # b3
            int.mat[,3]/sd1^3-int.mat[,1]/sd1)                    # sd1
  A = A*prop

  # N*4  integrate of derivative of log(P(s=1)) over y and d within (-Inf,y1)U(y2,Inf)
  B = cbind(-1*(r1*(r5d[,"y11"]-r5d[,"y01"]+r5d[,"y31"]-r5d[,"y21"])+(1-r1)*(r5d[,"y12"]-r5d[,"y02"]+r5d[,"y32"]-r5d[,"y22"])),   # t1
            r2*(r5p[,"y11"]-r5p[,"y01"]+r5p[,"y31"]-r5p[,"y21"])-r2*(r5p[,"y12"]-r5p[,"y02"]+r5p[,"y32"]-r5p[,"y22"]),        # t2
            -1*(r1*(r5d[,"y11"]-r5d[,"y01"]+r5d[,"y31"]-r5d[,"y21"])),                            # b3
            r1*(r6[,"y11"]-r6[,"y01"]+r6[,"y31"]-r6[,"y21"])+(1-r1)*(r6[,"y12"]-r6[,"y02"]+r6[,"y32"]-r6[,"y22"]))            # sd1
  B = B*prop

  D = A-B
  f = cbind(1,GE)
  nf = ncol(f)
  # change from N*4 to N*(2+2*nf)
  Df = cbind(matrix(D[,2],N,nf)*f,    # g0,g1,g2,...
             matrix(D[,1],N,nf)*f,    # b0,b1,b2,...
             D[,3],                   # b3
             D[,4])                   # sd1
  fisher.0a = t(Df) %*% Df
  rsum = colSums(Df)
  fisher.0b = rsum %*% t(rsum)
  fisher.0 = fisher.0b - fisher.0a;

  # N*10  integrate of (derivative of log(P(y,d)))*(derivative of log(P(y,d))) over y and d within (-Inf,y1)U(y2,Inf)
  AA = matrix(NA,N,10)
  AA[,1] = int.mat[,3]/sd1^4;                                         # t1*t1
  AA[,2] = (int.mat.1[,2]*r2+int.mat.0[,2]*(-1)*r2)/sd1^2;            # t1*t2
  AA[,3] = int.mat.1[,3]/sd1^4*r1;                                    # t1*b3
  AA[,4] = int.mat[,4]/sd1^5-int.mat[,2]/sd1^3;                       # t1*sd1
  # AA[,5] = r2^2*(1/r1*int.mat.1[,1]+1/(1-r1)*int.mat.0[,1]);          # t2*t2
  AA[,5] = r2^2*(1/r1*int.mat.1[,1]+1/r1s*int.mat.0[,1]);          # t2*t2
  AA[,6] = int.mat.1[,2]*r2/sd1^2;                                    # t2*b3
  AA[,7] = (int.mat.1[,3]*r2+int.mat.0[,3]*(-1)*r2)/sd1^3-(int.mat.1[,1]*r2+int.mat.0[,1]*(-1)*r2)/sd1;  # t2*sd1
  AA[,8] = AA[,3]                                                     # b3*b3
  AA[,9] = (int.mat.1[,4]/sd1^5-int.mat.1[,2]/sd1^3)*r1;              # b3*sd1
  AA[,10] = int.mat[,5]/sd1^6-2*int.mat[,3]/sd1^4+int.mat[,1]/sd1^2   # sd1*sd1
  AA = AA*prop;

  fisher = matrix(0,2+2*nf,2+2*nf)
  for(i in 1:N){
    ge = GE[i,];
    a = c(1,ge)
    TR = rbind(cbind(0,a,0,0),
               cbind(a,0,0,0),
               c(rep(0,2),1,0),
               c(rep(0,2),0,1))
    fisher.1 = B[i,]%*%t(B[i,])-A[i,]%*%t(B[i,])-B[i,]%*%t(A[i,])
    fisher.2 = matrix(NA,4,4)
    fisher.2[lower.tri(fisher.2,diag=T)]=AA[i,];
    fisher.2[upper.tri(fisher.2)] = t(fisher.2)[upper.tri(fisher.2)]
    fisher = fisher + TR%*%(fisher.1+fisher.2)%*%t(TR)
  }
  fisher = fisher+fisher.0
  colnames(fisher) = rownames(fisher) = c("g0",paste("g1",1:n.GE,sep="-"),
                                          "b0",paste("b1",1:n.GE,sep="-"),
                                          "b3","sd1")
  return(fisher)
}



