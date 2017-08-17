library(dlm)

sample.dlm <- array(NA, dim=c(N,nLatent,T+1));

dlm.model<-dlm(m0=m0,C0=C0,FF=F,V=V,W=W,GG=G)

start.time <- Sys.time()

dlmFilt<-dlmFilter(t(ct),dlm.model,debug=FALSE,simplify=TRUE)

for (i in 1:N) {
    sample<-dlmBSample(dlmFilt)
    sample.dlm[i,,] <- t(sample)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

cmprhist <- function(x,t,nclass) {
    h1 <- hist(eval(parse(text=sprintf('sample.%d.%d',x,t))),nclass = nclass)
    h2 <- hist(sample.dlm[,x,t],nclass = nclass)
    plot( h1, col=rgb(0,0,1,1/4), ylim=c(0,30500))  # first histogram
    plot( h2, col=rgb(1,0,0,1/4), ylim=c(0,30500), add=T)  # second
}

## compare filtering variances
v<-dlmSvd2var(dlmFilt$U.C, dlmFilt$D.C)
for (i in 1:(dim(sample.dlm)[3])) {
    if (max(abs(v[[i]]-eval(parse(text=sprintf('cov.%d', i)))))>1e-10) {
        stop(sprintf('%dth variance comparison fails', i));
    }
}
print("filtering variance comparison passed")

## compare filtering means
if (max(abs(fltmean-dlmFilt$m)) > 1e-10) {
    stop(sprintf('the filtering mean comparison fails'));
}
print("the filtering mean comparison passed")

## compare samples marginal distribution
for (i in 1:nLatent) {
    for (j in 1:(T+1)) {
        s1 = eval(parse(text=sprintf('sample.%d.%d', i, j)))
        s2 = sample.dlm[,i,j]
        if (t.test(s1,s2)$p.value < 0.05) {
            warning(sprintf('marginal mean t-test fails (%d,%d)', i,j));
        }
        if (var.test(s1,s2)$p.value < 0.05) {
            warning(sprintf('marginal variance t-test fails (%d,%d)', i,j));
        }
    }
}

