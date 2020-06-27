
## this script estimates sensitivity & specificity for methods that Jianjian has tried on tissue-trait associations

### load functions
source('~/Dropbox/Desktop/Yudi/SeSp.R') # Xia implemented exact solution for 3 methods
source('~/Dropbox/Desktop/Yudi/sens-spec-xia.r') # Yudi's ML estimation

### exact solution for 3 methods
SeSp <- function(count) {
	# count: data.frame with 3 columns of '+' and '-' for t1, t2, t3, and an n column for counts
	N <- sum(count$n)
	p1 <- sum(count$n[count$t1 == '+'])/N
	p2 <- sum(count$n[count$t2 == '+'])/N
	p3 <- sum(count$n[count$t3 == '+'])/N
	p12 <- sum(count$n[count$t1 == '+' & count$t2 == '+'])/N
	p23 <- sum(count$n[count$t2 == '+' & count$t3 == '+'])/N
	p13 <- sum(count$n[count$t1 == '+' & count$t3 == '+'])/N
	p123 <- sum(count$n[count$t1 == '+' & count$t2 == '+' & count$t3 == '+'])/N
	V2 <- (p123 - p12*p3 - p13*p2 - p23*p1 + 2*p1*p2*p3)**2/((p12 - p1*p2)*(p13 - p1*p3)*(p23 - p2*p3))
	#rho <- .5 - sqrt(.25 - 1/(4 + V**2)) # or minus
	rho <- .5 - sqrt(.25 - 1/(4 + V2))
	#if (rho > 1 | rho < 0) rho <- .5 - sqrt(.25 - 1/(4 + V**2))
	C1 <- (p12 - p1*p2)*(p13 - p1*p3)/(p23 - p2*p3)
	C2 <- (p12 - p2*p1)*(p23 - p2*p3)/(p13 - p1*p3)
	C3 <- (p13 - p3*p1)*(p23 - p3*p2)/(p12 - p1*p2)
	sens1 <- tp1 <- p1 + sqrt(C1)*sqrt((1 - rho)/rho)
	fp1 <- p1 - sqrt(C1)*sqrt(rho/(1 - rho))
	spec1 <- 1 - fp1
	sens2 <- tp2 <- p2 + sqrt(C2)*sqrt((1 - rho)/rho)
	fp2 <- p2 - sqrt(C2)*sqrt(rho/(1 - rho))
	spec2 <- 1 - fp2
	sens3 <- tp3 <- p3 + sqrt(C3)*sqrt((1 - rho)/rho)
	fp3 <- p3 - sqrt(C3)*sqrt(rho/(1 - rho))
	spec3 <- 1 - fp3
	ee <- c(sens1, sens2, sens3, spec1, spec2, spec3, rho)
	ee[ee > 1] <- 1
	return(ee)
}

## yudi's function

mixprob = function(sens,spec,p)  # 
{
	ntest = length(sens)
	a=list()
	for (i in 1:ntest) a[[i]]=0:1
	ytab = expand.grid(a)
	cname = apply(ytab,1,paste0, collapse='')
	##
	prob1 = t(sens**t(ytab)*(1-sens)**t(1-ytab))
	prob1 = apply(prob1,1,prod)
	#prob1 = prob1/sum(prob1)
	##
	prob0 = t((1-spec)**t(ytab)*spec**t(1-ytab))
	prob0 = apply(prob0,1,prod)
	#prob0 = prob0/sum(prob0)
	mixprob = p*prob1 + (1-p)*prob0 
	names(mixprob) = cname
	return(mixprob)
}

##
## ny = summary counts vector, must be named according 0-1 test result
## eg 010 = results from 3 tests are 0,1,0 respectively
## NOTE: potential complications:
##   names(ny) maybe incomplete and has arbitrary order
##
## return MINUS loglik for minimization!!! 
loglik = function(param, ntest, ny) #, p = .1) 
{
	if (min(param)<=0.0001| max(param)>=1-0.0001) {loglik = exp(30);
		return(loglik)}
	
	sens= param[1:ntest]
	spec = param[(ntest+1):(2*ntest)]
	p = param[2*ntest+1]
	prob = mixprob(sens,spec,p)  ## this has complete list of names
	##   
	obs.comb = intersect(names(ny), names(prob))
	obs.ny = ny[obs.comb]
	model.prob = prob[obs.comb]
	if (min(prob)>0.0001 & max(prob)<1-0.0001) loglik = -sum(obs.ny*log(model.prob))
	if (min(prob)<=0.0001| max(prob)>=1-0.0001) loglik = exp(30)
	
	return(loglik)
}

### load Jianjian's results
#load('~/Dropbox/Desktop/Jianjian/three_matrix2Xia.Rdata')
load('~/Dropbox/Desktop/Jianjian/four_methods_pmatrix2Xia.Rdata')
dd <- data.frame(as.numeric(as.matrix(eQTL)), as.numeric(as.matrix(LDSC)), as.numeric(as.matrix(LDSC_100K)), as.numeric(as.matrix(ROLYPOLY)), as.numeric(as.matrix(LM)), as.numeric(as.matrix(GLM)))

pv <- data.frame(eqtl = as.numeric(as.matrix(eQTL)), ldsc = as.numeric(as.matrix(LDSC)), rolypoly = as.numeric(as.matrix(ROLYPOLY)), jianjian = as.numeric(LM)) #, jianjian2 = as.numeric(GLM))
pv <- na.omit(pv)

pthres <- c(10**(-(16:1)), seq(.1, 1, length = 20))
pthres[16] <- .05

# 1e-5, 1e-4, 5e-4
# 1e-5, 5e-4, 1e-3
# 5e-4, 5e-3, 5e-2



thres.eqtl <- rep(5e-4, nrow(pv))
thres.ldsc <- rep(5e-4, nrow(pv))
thres.jianjian <- rep(5e-4, nrow(pv))
thres.rolypoly <- rep(5e-4, nrow(pv))
thres <- data.frame(eqtl = thres.eqtl, ldsc = thres.ldsc, rolypoly = thres.rolypoly, jianjian = thres.jianjian)

# res <- res1 <- matrix(NA, length(pthres), 7)
# count <- data.frame(t1 = c('+', '+', '+', '+', '-', '-', '-', '-'), 
#  		            t2 = c('+', '+', '-', '-', '+', '+', '-', '-'),
#  					t3 = c('+', '-', '+', '-', '+', '-', '+', '-'), 
#  					n = c(207, 27, 75, 31, 50, 29, 85, 162))
# x  <- (pv < .1)*1
# xtab <- table(paste0(x[,1], x[,2], x[,3]))
# xtab[1:8] <- 0

#require(svMisc)
#for (i in 1:length(pthres)) {
	xx  <- (pv < thres)*1
	xxtab <- table(paste0(xx[,1], xx[,2], xx[,3], xx[,4]))
	#xtab[names(xxtab)] <- xxtab
	#count$n <- rev(xtab)
	#res[i,] <- SeSp(count)
	ny.ran = xxtab
	#ny.ran[1] <- 0 ## try no 000
	param = c(rep(.5, 8), .5)
	res1 = optim(param, loglik, ntest = 4, ny = ny.ran)$par ## Nelder-Mead
#	progress(i/length(pthres)*100)
#}

res.boot <- matrix(NA, 99, 9)
require(svMisc)
for (j in 1:99) {
	#cat(j, '\n')
	pv2 <- pv[sample(1:nrow(pv), nrow(pv), replace = TRUE),]
	#for (i in 1:length(pthres)) {
		#i = 13
		xx  <- (pv2 < thres)*1
		xxtab <- table(paste0(xx[,1], xx[,2], xx[,3], xx[,4]))
		#xtab[names(xxtab)] <- xxtab
		#count$n <- rev(xtab)
		ny.ran = xxtab
		#ny.ran[1] <- 0 ## try no 000
		param = c(rep(.5, 8), .5)
		res.boot[j,] = optim(param, loglik, ntest = 4, ny = ny.ran)$par ## Nelder-Mead
	progress(j/99*100)
	#}
	#cat('\n')
}

rb <- res.boot
rb <- na.omit(rb)

idx <- which(rb[,2] < 1 - rb[,5])
tmp <- 1 - rb[idx,2] 
rb[idx,2] <- 1 - rb[idx,5]
rb[idx,5] <- tmp
tmp <- 1 - rb[idx,3] 
rb[idx,3] <- 1 - rb[idx,6]
rb[idx,6] <- tmp
#idx <- which(rb[,2] < .5 | rb[,3] < .5)
#rb <- rb[-idx,]

require(RColorBrewer)
co <- brewer.pal(4, "Set1")
plot(0:1, ann = FALSE, axes = FALSE, type = 'n', xlim = c(0, 1), ylim = c(0, 1))
axis(1, at = seq(0, 1, .2), labels = seq(1, 0, -.2))
mtext('Specificity', 1, 3)
axis(2, at = seq(0, 1, .2), labels = seq(0, 1, .2))
mtext('Sensitivity', 2, 3)
abline(0, 1, lty = 3)
points(1 - res1[5], res1[1], col = paste0(co[1]), pch = 16, cex = 1.5)
points(1 - res1[6], res1[2], col = paste0(co[2]), pch = 16, cex = 1.5)
points(1 - res1[7], res1[3], col = paste0(co[3]), pch = 16, cex = 1.5)
points(1 - res1[8], res1[4], col = paste0(co[4]), pch = 16, cex = 1.5)
points(1 - rb[,5], rb[,1], col = paste0(co[1], '4D'), pch = 16)
points(1 - rb[,6], rb[,2], col = paste0(co[2], '4D'), pch = 16)
points(1 - rb[,7], rb[,3], col = paste0(co[3], '4D'), pch = 16)
points(1 - rb[,8], rb[,4], col = paste0(co[4], '4D'), pch = 16)
legend(1 - .2, .2, c('eQTL', 'LDSC', 'RolyPoly', 'Zhijian'), pch = c(16, 16, 16, 16), col = co, bty = 'n')






## big simu data processing

# .5, .5, .5, .5 etc. replot

for (th1 in c(5e-1, 5e-2, 5e-3, 5e-4)) {
	cat(th1, '\n')
	for (th2 in c(5e-1, 5e-2, 5e-3, 5e-4)) {
		for (th3 in c(5e-1, 5e-2, 5e-3, 5e-4)) {
			for (th4 in c(5e-1, 5e-2, 5e-3, 5e-4)) {
				pdf(paste0('~/Dropbox/Desktop/Jianjian/4tests_bootstrap/replot_rb_', th1, '_', th2, '_', th3, '_', th4, '.pdf'), width = 4.2, height = 3)
				d <- readRDS(paste0('~/Dropbox/Desktop/Jianjian/4tests_bootstrap/rb_', th1, '_', th2, '_', th3, '_', th4, '.rds'))
				data <- data.frame(TPR = as.numeric(d[,1:4]), FPR = 1 - as.numeric(d[,5:8]), method = rep(c('eQTL', 'LDSC', 'RolyPoly', 'ChiSpec'), each = 100))
				require(ggplot2)
				require(RColorBrewer)
				co <- brewer.pal(9, "Set1")[c(5,3,2,1)]
				p <- ggplot(data, aes(x = FPR, y = TPR)) + 
				  stat_density2d(geom = "tile", aes(fill = method, alpha = ..density..), contour = FALSE, n = 100, h = c(.36, .36)) + 
				  scale_alpha(range = c(0, 1)) + 
				  scale_fill_manual(values=c("eQTL"=co[1], "LDSC"=co[2], "RolyPoly"=co[3], "ChiSpec"=co[4])) + #geom_point() +
				  theme_minimal() +
				  xlim(0, 1) + ylim(0, 1) 
				print(p)
				dev.off()
			}
		}
	}
}


### compute mean sens, spec for each test/threshold

est <- matrix(NA, 4, 9) # 4 thresholds, 4 tests, 9 parameters
dimnames(est) <- list(c(5e-1, 5e-2, 5e-3, 5e-4), c('sens.eQTL', 'sens.LDSC', 'sens.RolyPoly', 'sens.Zhijian', 'spec.eQTL', 'spec.LDSC', 'spec.RolyPoly', 'spec.Zhijian', 'prev'))
se <- est

### load the data

data <- c()	
for (th1 in c(5e-1, 5e-2, 5e-3, 5e-4)) {
	cat(th1, '\n')
	for (th2 in c(5e-1, 5e-2, 5e-3, 5e-4)) {
		for (th3 in c(5e-1, 5e-2, 5e-3, 5e-4)) {
			for (th4 in c(5e-1, 5e-2, 5e-3, 5e-4)) {
				d <- readRDS(paste0('~/Dropbox/Desktop/Jianjian/4tests_bootstrap/rb_', th1, '_', th2, '_', th3, '_', th4, '.rds'))
				df <- data.frame(th1 = th1, th2 = th2, th3 = th3, th4 = th4, d)
				data <- rbind(data, df)
			}
		}
	}
	cat('\n')
}

### load the data while computing est + se

data.mean <- data.var <- c()
for (th1 in c(5e-1, 5e-2, 5e-3, 5e-4)) {
	cat(th1, '\n')
	for (th2 in c(5e-1, 5e-2, 5e-3, 5e-4)) {
		for (th3 in c(5e-1, 5e-2, 5e-3, 5e-4)) {
			for (th4 in c(5e-1, 5e-2, 5e-3, 5e-4)) {
				d <- readRDS(paste0('~/Dropbox/Desktop/Jianjian/4tests_bootstrap/rb_', th1, '_', th2, '_', th3, '_', th4, '.rds'))
				bad <- which(d[,1] < 1 - d[,5] | d[,2] < 1 - d[,6] | d[,3] < 1 - d[,7] | d[,4] < 1 - d[,8]) # remove dots below diagonal
				d <- d[-bad,]
				if (is.null(nrow(d))) d <- rbind(d, d)
				means <- colMeans(d)
				vars <- colVars(d)/nrow(d) # removed dots below diagonal
				df <- data.frame(th1 = th1, th2 = th2, th3 = th3, t(means))
				data.mean <- rbind(data.mean, df)
				df <- data.frame(th1 = th1, th2 = th2, th3 = th3, t(vars))
				data.var <- rbind(data.var, df)
			}
		}
	}
	cat('\n')
}

### meta the means

data.mean <- na.omit(data.mean)
data.var <- na.omit(data.var)
bad <- which(data.var[,4] == 0)
data.mean <- data.mean[-bad,]
data.var <- data.var[-bad,]
for (th in c(5e-1, 5e-2, 5e-3, 5e-4)) {
	for (i in 1:9) {
		mt <- meta(data.mean[data.mean$th1 == th,i + 3], sqrt(data.var[data.var$th1 == th,i + 3]))
		est[as.character(th),i] <- mt$beta
		se[as.character(th),i] <- mt$se
	}
}

### for 0.05, 0.01
for (th in c(5e-1, 5e-2, 5e-3, 5e-4)) {
	for (i in 1:9) {
		mt <- meta(data.mean[data.mean$th1 == th,i + 3], sqrt(data.var[data.var$th1 == th,i + 3]))
		est[as.character(th),i] <- mt$beta
		se[as.character(th),i] <- mt$se
	}
}

### ROC

require(RColorBrewer)
co <- brewer.pal(4, "Set1")

plot(0:1, ann = FALSE, axes = FALSE, type = 'n', xlim = c(0, 1), ylim = c(0, 1))
axis(1, at = seq(0, 1, .2), labels = seq(1, 0, -.2))
mtext('Specificity', 1, 3)
axis(2, at = seq(0, 1, .2), labels = seq(0, 1, .2))
mtext('Sensitivity', 2, 3)
abline(0, 1, lty = 3)
points(1 - est[,5], est[,1], col = paste0(co[1]), pch = 16)
points(1 - est[,6], est[,2], col = paste0(co[2]), pch = 16)
points(1 - est[,7], est[,3], col = paste0(co[3]), pch = 16)
points(1 - est[,8], est[,4], col = paste0(co[4]), pch = 16)
points(1 - data.mean[,5 + 3], data.mean[,1 + 3], col = paste0(co[1], '4D'), pch = 16, cex = .36)
points(1 - data.mean[,6 + 3], data.mean[,2 + 3], col = paste0(co[2], '4D'), pch = 16, cex = .36)
points(1 - data.mean[,7 + 3], data.mean[,3 + 3], col = paste0(co[3], '4D'), pch = 16, cex = .36)
points(1 - data.mean[,8 + 3], data.mean[,4 + 3], col = paste0(co[4], '4D'), pch = 16, cex = .36)
legend(1 - .2, .2, c('eQTL', 'LDSC', 'RolyPoly', 'Zhijian'), pch = c(16, 16, 16, 16), col = co, bty = 'n')


### TDR vs #discovery
prev = data.mean[,7 + 3]
tp1 = prev*data.mean[,1 + 3]
tp2 = prev*data.mean[,2 + 3]
tp3 = prev*data.mean[,3 + 3]
tn1 = (1 - prev)*data.mean[,4 + 3]
tn2 = (1 - prev)*data.mean[,5 + 3]
tn3 = (1 - prev)*data.mean[,6 + 3]
fn1 = prev - tp1
fn2 = prev - tp2
fn3 = prev - tp3
fp1 = (1 - prev) - tn1
fp2 = (1 - prev) - tn2
fp3 = (1 - prev) - tn3

tdr1 = tp1/(tp1 + fp1)
tdr2 = tp2/(tp2 + fp2)
tdr3 = tp3/(tp3 + fp3)

disc1 = tp1 + fp1
disc2 = tp2 + fp2
disc3 = tp3 + fp3

co <- ggcolor(10)
plot(0:1, ann = FALSE, axes = FALSE, type = 'n', xlim = c(0, 1), ylim = c(0, 1))
axis(1, at = seq(0, 1, .2), labels = seq(0, 1, .2))
mtext('# Positive', 1, 3)
axis(2, at = seq(0, 1, .2), labels = seq(0, 1, .2))
mtext('TDR', 2, 3)
abline(0, 1, lty = 3)
points(disc1, tdr1, col = paste0(co[1], '4D'), pch = 16, cex = .36)
points(disc2, tdr2, col = paste0(co[4], '4D'), pch = 16, cex = .36)
points(disc3, tdr3, col = paste0(co[7], '4D'), pch = 16, cex = .36)
legend(1 - .2, .2, c('eQTL', 'LDSC', 'Zhijian'), pch = c(16, 16, 16), col = co[c(1, 4, 7)], bty = 'n')

prev = est[,7]
tp1 = prev*est[,1]
tp2 = prev*est[,2]
tp3 = prev*est[,3]
tn1 = (1 - prev)*est[,4]
tn2 = (1 - prev)*est[,5]
tn3 = (1 - prev)*est[,6]
fn1 = prev - tp1
fn2 = prev - tp2
fn3 = prev - tp3
fp1 = (1 - prev) - tn1
fp2 = (1 - prev) - tn2
fp3 = (1 - prev) - tn3

tdr1 = tp1/(tp1 + fp1)
tdr2 = tp2/(tp2 + fp2)
tdr3 = tp3/(tp3 + fp3)

disc1 = tp1 + fp1
disc2 = tp2 + fp2
disc3 = tp3 + fp3

points(disc1, tdr1, col = paste0(co[1]), pch = 16)
points(disc2, tdr2, col = paste0(co[4]), pch = 16)
points(disc3, tdr3, col = paste0(co[7]), pch = 16)
















