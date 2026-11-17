suppressMessages({library(ape);library(phylolm)})
tr<-read.tree('phaseC/f311_GBM_t3_chronogram.nwk');d<-read.csv('phaseC/len311.csv');rownames(d)<-d$tip
d<-d[d$tip %in% tr$tip.label,];tr<-keep.tip(tr,d$tip)
for(m in c('len~nfun','len~nfun+para+myco','len~nfun*carn+para+myco')){f<-phylolm(as.formula(m),d,tr,model='lambda');print(m);print(summary(f)$coefficients);cat('lambda',f$optpar,'n',nrow(d),'\n')}
