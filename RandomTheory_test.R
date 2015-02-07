library(MASS)
library(splines)
library(Matrix)
library(lattice)
library(expm)
library(foreach)
library(doSNOW)
library(lme4)


num.repli = 100
n = 50
J = 10
N = n*J
pho.x = 0.4

gg = 1
aa = 1
bb = 0
ss = 0.4
pho.x = 0.4


gamma0Fun = function(t){2*t^2 - 1}
gammaFun = function(t){
	gamma1 = gg*(t > 0.4)
	gamma2 = -gg*cos(2*pi*t)
	gamma3 = gg*((2 - 3*t)^2/2 - 1)
	gamma4 = gg*sin(2*pi*t)
	return(c(gamma1, gamma2, gamma3, gamma4))
			}

alphaFun = function(t){
	alpha1 = aa*(cos(2*pi*t)+bb)
	alpha2 = aa*(sin(2*pi*t)+bb)
	alpha3 = aa*(exp(-2*t/(t+1))+bb)
	return(c(alpha1,alpha2,alpha3))
}

q = 3			
sigma.random = matrix(c(6,3.8,0.6, 3.8, 5, 1, 0.6, 1, 4), q, q)*ss

p.all = 100
index.nz.r = c(2, 50, 90)
pho.random = 0.4

nknots = 1
num.degree = 3

## Covariance matrix for generating X

sigma.x = 1
cov.x = diag(1, p.all+1)
cov.x = sigma.x^2 * pho.x^abs(row(cov.x) - col(cov.x))
ev.x = eigen(cov.x, symmetric = TRUE)
   if (!all(ev.x$values >= -sqrt(.Machine$double.eps) * abs(ev.x$values[1]))) {
            warning("sigma is numerically not positive definite")
        }
cov.x.root = ev.x$vectors %*% diag(sqrt(ev.x$values), length(ev.x$values)) %*% t(ev.x$vectors)


ncl = 30
TYPE = "SOCK"
cluster = makeCluster(ncl, type = TYPE)
registerDoSNOW(cluster)

Sigma = Sigma.raw = matrix(0, q*J, q*J)

Row = sigma.random
for(k in 1:(J-1)){
	Row = cbind(Row, sigma.random)
			}
Sigma.raw = Row
for(k in 1:(J-1)){
	Sigma.raw = rbind(Sigma.raw, Row)
			}
## Function to generate covariance matrix for random effect for each subject
RanCov = function(t){
	for(i in 1:(q*J)){
		for(j in 1:(q*J)){
			diff = abs(t[ceiling(i/q)] - t[ceiling(j/q)])
			Sigma[i, j] = Sigma.raw[i, j]*(pho.random^diff)
						}
					}
	return(Sigma)	
}

## Function used to generate basis for beta1
Beta1BsGen = function(A_comb){
	return(A_comb[-1]*A_comb[1])
					}

J_simu = rep(J, n)
id = c(t(matrix(rep(1:n, J), nrow = n)))

clusterExport(cluster, c("Sigma.raw","Sigma","J_simu","id"))


RESULT_ALL = foreach(repli = 1:num.repli, .combine = "rbind", .packages = c("lme4","foreach","MASS","splines")) %dopar% {

	set.seed(2014+repli)

tstar.x = foreach(i=1:n, .combine = "rbind") %dopar% matrix(rnorm(J*(p.all+1)), nrow = J, ncol = p.all+1) %*% cov.x.root


time.all.onelist = pnorm(tstar.x[,1])
	covar.all = tstar.x[,-1]
	
	randomcoeff.all = foreach(i = 1:n, .combine = "rbind",.packages = "MASS") %dopar%{
		if (i == 1) start = 0 else {start = sum(J_simu[1:(i-1)])}
		end = sum(J_simu[1:i])
		idx_i = (start+1):end
		mvrnorm(1, rep(0,q*J), RanCov(time.all.onelist[idx_i]))
	}
	
pred.all = covar.all[,index.nz.r]	

effect.all = foreach(i = 1:n, .combine = c, .packages = "foreach") %dopar%{
		
		if (i == 1) start = 0 else {start = sum(J_simu[1:(i-1)])}
		end = sum(J_simu[1:i])
		idx_i = (start + 1):end
		
		B = matrix(randomcoeff.all[i, ], nrow = J, ncol = q, byrow = TRUE)
		A = t(sapply(time.all.onelist[idx_i], alphaFun))
		rancoef_i = foreach(b = 1:J, .combine = "rbind") %dopar% A[b, ]*B[b,]
		
		pred_i = pred.all[idx_i,]
		
		Y_i = foreach(a = 1:J, .combine = c) %dopar% sum(pred_i[a, ]*rancoef_i[a,])
		}
	Y.all.onelist = effect.all + rnorm(N)
	
bs.beta0 = cbind(rep(1, length(time.all.onelist)), bs(time.all.onelist, degree = num.degree, df = 1 + nknots + num.degree))

	
random_result = foreach(k = 1:p.all, .combine = "rbind",.packages = c("lme4")) %dopar% {
	bs.random = t(apply(cbind(covar.all[,k], bs.beta0), 1, Beta1BsGen))
	X = cbind(bs.beta0, bs.random)	
	fs.lmer = lmer(Y.all.onelist ~ (0 + X[,7]|id) + (0 + X[,8]|id) + (0 + X[,9]|id) + (0 + X[,10]|id) + (0 + X[,11]|id) + (0 + X[,12]|id),data = as.data.frame(X))
	
	pss = sum(fitted(fs.lmer)^2)
	rss = sum(residuals(fs.lmer)^2)

	Dk = diag(VarCorr(fs.lmer))
	
	theory_cri = foreach(i =1:n, .combine = "rbind")%dopar%{
		
	if (i == 1) start = 0 else {start = sum(J_simu[1:(i-1)])}
	end = sum(J_simu[1:i])
	idx_i = (start + 1):end
	
	Uk_i = X[idx_i,7:12]
	Ak_i = Uk_i%*%Dk%*%t(Uk_i) 
	Vk_i = Ak_i + attr(VarCorr(fs.lmer), "sc")^2*diag(1, J)
	c1 = sum(diag((Ak_i)))
	c2 = sum(diag((Ak_i %*% solve(Vk_i) %*% Ak_i)))
	c(c1, c2)
	}
	tc = apply(theory_cri, 2, mean)
		
	c(pss, rss, tc)
}
rank_pss = p.all + 1 - rank(random_result[,1])
rank_rss = rank(random_result[,2])
rank_tc1 = p.all + 1 - rank(random_result[,3])
rank_tc2 = p.all + 1 - rank(random_result[,4])

names(rank_pss) = c()
names(rank_rss) = c()
names(rank_tc1) = c()
names(rank_tc2) = c()

c(rank_pss[index.nz.r],rank_rss[index.nz.r], rank_tc1[index.nz.r], rank_tc2[index.nz.r])

}



theory_V2 = foreach(i =1:n)%dopar%{
		
	if (i == 1) start = 0 else {start = sum(J_simu[1:(i-1)])}
	end = sum(J_simu[1:i])
	idx_i = (start + 1):end
	
	Uk_i = X[idx_i,7:12]
	Ak_i = Uk_i%*%Dk%*%t(Uk_i) 
	Vk_i = Ak_i + attr(VarCorr(fs.lmer), "sc")^2*diag(1, J)
	}
a1 = theory_V1[[2]]
a2 = theory_V2[[2]]
round(a1[1:5,1:5], 3);round(a2[1:5,1:5], 3);
round(theory_V1[[1]],4)