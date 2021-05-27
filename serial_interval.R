library(sn)
library(readr)
library(data.table)
library(fitdistrplus)

setwd('C:/Users/Sony/Desktop/VOC617')

casePair_VOC <- fread('householdPair_b.1.617.2.csv')
casePair_preLD <- fread('householdPair_preLD.csv')

# sample serial interval
t=sort(casePair_preLD[,unique(ONSET_TO_ISOLATION_INFECTOR)])

listSerial = lapply(1:(max(t)+1), function(i){
  casePair_preLD[ONSET_TO_ISOLATION_INFECTOR==i-1, SERIAL_INT]
})

sampleSerialInt = data.table(ITER = rep(1:1000, each = casePair_VOC[,.N]),
                             ONSET_TO_ISOLATION = rep(casePair_VOC[,ONSET_TO_ISOLATION_INFECTOR], times = 1000))
sampleSerialInt[, SERIAL_INTERVAL := sample(listSerial[[(ONSET_TO_ISOLATION+1)]], 1, replace = T), 
                by = seq_len(nrow(sampleSerialInt))]

# test fit parameteric distribution 
dat = sampleSerialInt[ITER == 1, SERIAL_INTERVAL]
skew_mod <-  selm(dat ~ 1) # selm is "skew-elliptic lm"
summary(skew_mod)

hist(dat,prob=TRUE,nclass="scott") # "scott" is from MASS
plot(function(x) dsn(x, dp=skew_mod@param$dp), from=-3, to=14, col="red", add=TRUE)

# fit parametric distribution for each iteration's samples
set.seed(123)
fitSerialInt = lapply(1:1000, function(i){
  
  x = sampleSerialInt[ITER == i, SERIAL_INTERVAL]
  skew_mod =  selm(x ~ 1) # selm is "skew-elliptic lm"
  fit = dsn(-3:14, dp=skew_mod@param$dp)
  fit = fit/sum(fit)
  
})

distSerialInt = data.table(ITER = rep(1:1000, each = length(-3:14)), 
                           SERIAL_INTERVAL = rep(-3:14, times = 1000),
                           PDF = unlist(fitSerialInt))
distSerialInt[,CDF:=cumsum(PDF), by=.(ITER)]

distSerialInt[sampleSerialInt[,.N, by=.(ITER,SERIAL_INTERVAL)],
              PDF_NON_FIT:=i.N,on=c(ITER='ITER',SERIAL_INTERVAL='SERIAL_INTERVAL')]
distSerialInt[is.na(PDF_NON_FIT),PDF_NON_FIT:=0]
distSerialInt[,PDF_NON_FIT:=PDF_NON_FIT/casePair_VOC[,.N]]
distSerialInt[,CDF_NON_FIT:=cumsum(PDF_NON_FIT),by=.(ITER)]

# mean and 95%CI of mean serial interval for all iterations
distSerialInt[,MEAN:=PDF*SERIAL_INTERVAL]
summary(distSerialInt[,sum(MEAN), by=.(ITER)]$V1)
quantile(distSerialInt[,sum(MEAN), by=.(ITER)]$V1, probs = c(0.025,0.975))

# tabulate median serial interval for all iterations
mean(distSerialInt[CDF>=0.5, min(SERIAL_INTERVAL), by=.(ITER)]$V1)
table(distSerialInt[CDF>=0.5, min(SERIAL_INTERVAL), by=.(ITER)]$V1)
quantile(distSerialInt[CDF>=0.5, min(SERIAL_INTERVAL), by=.(ITER)]$V1, probs = c(0.025,0.975))

# tabulate mode serial interval for all iterations
t = -3:14
t[distSerialInt[,which.max(PDF), by=.(ITER)]$V1]
mean(t[distSerialInt[,which.max(PDF), by=.(ITER)]$V1])
table(t[distSerialInt[,which.max(PDF), by=.(ITER)]$V1])
quantile(t[distSerialInt[,which.max(PDF), by=.(ITER)]$V1], probs = c(0.025,0.975))

distSerialInt = distSerialInt[,.(mean(PDF), quantile(PDF, probs = 0.025), quantile(PDF, probs = 0.975),
                                 mean(CDF), quantile(CDF, probs = 0.025), quantile(CDF, probs = 0.975),
                                 mean(PDF_NON_FIT), quantile(PDF_NON_FIT, probs = 0.025), quantile(PDF_NON_FIT, probs = 0.975),
                                 mean(CDF_NON_FIT), quantile(CDF_NON_FIT, probs = 0.025), quantile(CDF_NON_FIT, probs = 0.975)), 
                              by=.(SERIAL_INTERVAL)]
setnames(distSerialInt, old=c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12'), 
         new=c('PMF_MEAN', 'PMF_LOWER_CI', 'PMF_UPPER_CI', 'CMF_MEAN', 'CMF_LOWER_CI', 'CMF_UPPER_CI',
               'PMF_NF_MEAN', 'PMF_NF_LOWER_CI', 'PMF_NF_UPPER_CI', 'CMF_NF_MEAN', 'CMF_NF_LOWER_CI', 'CMF_NF_UPPER_CI'))
