library(sn)
library(readr)
library(data.table)
library(fitdistrplus)
set.seed(321)

setwd('C:/Users/Sony/Desktop/VOC617/serial interval')

# load data
casePair.VOC <- fread('householdPair_b.1.617.2.csv') # household transmission pair infected with B.1.617.2 variant
casePair.preLD <- fread('householdPair_preLD.csv') # household transmission pair infected prior to partial lockdown

# prepare list of serial interval observed prior to partial lockdown stratified by duration from onset to isolation 
t=sort(casePair.preLD[,unique(ONSET_TO_ISOLATION_INFECTOR)])

listSerial = lapply(1:(max(t)+1), function(i){
  casePair.preLD[ONSET_TO_ISOLATION_INFECTOR==i-1, SERIAL_INT]
})

# for each duration of onset to isolation observed in recent B.1.617.2 household infector, 
# sample the equivalent serial interval observed prior to the partial lockdown and calculate the differences
sampleSerialInt = data.table(ITER = rep(1:1000, each = casePair.VOC[,.N]),
                             ONSET_TO_ISOLATION = rep(casePair.VOC[,ONSET_TO_ISOLATION_INFECTOR], times = 1000))
sampleSerialInt[, SERIAL_INTERVAL_PRE_LD := sample(listSerial[[(ONSET_TO_ISOLATION+1)]], 1, replace = T), 
                by = seq_len(nrow(sampleSerialInt))]


# test fit parameteric distribution 
dat = sampleSerialInt[ITER == 1, SERIAL_INTERVAL_PRE_LD]
skew_mod <-  selm(dat ~ 1) # selm is "skew-elliptic lm"
summary(skew_mod)

hist(dat,prob=TRUE,breaks = seq(-3,15,1))
plot(function(x) dsn(x, dp=skew_mod@param$dp), from=-3, to=15, col="red", add=TRUE)

# fit parametric distribution for each iteration's samples
# calculate pmf and cdf
fitSerialInt.preLD = lapply(1:1000, function(i){
  
  x = sampleSerialInt[ITER == i, SERIAL_INTERVAL_PRE_LD]
  skew_mod =  selm(x ~ 1) # selm is "skew-elliptic lm"
  fit = dsn(-3:15, dp=skew_mod@param$dp)
  fit = fit/sum(fit)
  
  return(list(skew_mod@param$dp, fit))
  
})

distSerialInt = data.table(ITER = rep(1:1000, each = length(-3:15)), 
                           SERIAL_INTERVAL = rep(-3:15, times = 1000),
                           PMF_PRE_LD = unlist(lapply(fitSerialInt.preLD, `[[`, 2)))
distSerialInt[,CDF_PRE_LD:=cumsum(PMF_PRE_LD), by=.(ITER)]



# mean and 95%CI of the mean of serial intervals sampled from pre lockdown cases for all iterations
distSerialInt[,MEAN_PRE_LD:=PMF_PRE_LD*SERIAL_INTERVAL]
mean_mean.preLD = mean(distSerialInt[,sum(MEAN_PRE_LD), by=.(ITER)]$V1)
ci_mean.preLD = quantile(distSerialInt[,sum(MEAN_PRE_LD), by=.(ITER)]$V1, probs = c(0.025,0.975))

# mean serial  interval of B.1.617.2 cases 
mean_VOC = mean(casePair.VOC$SERIAL_INT)

# differences in mean
difference.in.mean = mean_VOC - distSerialInt[,sum(MEAN_PRE_LD), by=.(ITER)]$V1
mean_difference.in.mean = mean(difference.in.mean)
ci_difference.in.mean = quantile(difference.in.mean, probs = c(0.025,0.975))

output_mean = c(mean_VOC, mean_mean.preLD, ci_mean.preLD, mean_difference.in.mean, ci_difference.in.mean)

# mean and 95%CI of the median of serial intervals sampled from pre lockdown cases for all iterations
mean_median.preLD = mean(distSerialInt[CDF_PRE_LD>=0.5, min(SERIAL_INTERVAL), by=.(ITER)]$V1)
table(distSerialInt[CDF_PRE_LD>=0.5, min(SERIAL_INTERVAL), by=.(ITER)]$V1)
ci_median.preLD = quantile(distSerialInt[CDF_PRE_LD>=0.5, min(SERIAL_INTERVAL), by=.(ITER)]$V1, probs = c(0.025,0.975))

# median serial  interval of B.1.617.2 cases 
median_VOC = median(casePair.VOC$SERIAL_INT)

# differences in median
difference.in.median = median_VOC - distSerialInt[CDF_PRE_LD>=0.5, min(SERIAL_INTERVAL), by=.(ITER)]$V1 
mean_difference.in.median = mean(difference.in.median)
ci_difference.in.median = quantile(difference.in.median, probs = c(0.025,0.975))

output_median = c(median_VOC, mean_median.preLD, ci_median.preLD, mean_difference.in.median, ci_difference.in.median)

# mean and 95%CI of the mode of serial intervals sampled from pre lockdown cases for all iterations
t = -3:15
t[distSerialInt[,which.max(PMF_PRE_LD), by=.(ITER)]$V1]
mean_mode.preLD = mean(t[distSerialInt[,which.max(PMF_PRE_LD), by=.(ITER)]$V1])
table(t[distSerialInt[,which.max(PMF_PRE_LD), by=.(ITER)]$V1])
ci_mode.preLD = quantile(t[distSerialInt[,which.max(PMF_PRE_LD), by=.(ITER)]$V1], probs = c(0.025,0.975))

# mode of the serial  interval distribution of B.1.617.2 cases 
mode_VOC = 2
  
# differences in mode
difference.in.mode = mode_VOC - t[distSerialInt[,which.max(PMF_PRE_LD), by=.(ITER)]$V1] 
mean_difference.in.mode = mean(difference.in.mode)
ci_difference.in.mode = quantile(difference.in.mode, probs = c(0.025,0.975))

output_mode = c(mode_VOC, mean_mode.preLD, ci_mode.preLD, mean_difference.in.mode, ci_difference.in.mode)


# prep distribution data for plots
distSerialInt = distSerialInt[,.(mean(PMF_PRE_LD), quantile(PMF_PRE_LD, probs = 0.025), quantile(PMF_PRE_LD, probs = 0.975),
                                 mean(CDF_PRE_LD), quantile(CDF_PRE_LD, probs = 0.025), quantile(CDF_PRE_LD, probs = 0.975)),
                              by=.(SERIAL_INTERVAL)]
setnames(distSerialInt, 
         c('SERIAL_INTERVAL', 'PMF_PRE_LD_MEAN', 'PMF_PRE_LD_LOWER_CI', 'PMF_PRE_LD_UPPER_CI', 
           'CDF_PRE_LD_MEAN', 'CDF_PRE_LD_LOWER_CI', 'CDF_PRE_LD_UPPER_CI'))


# prep summary table
outputTable = rbind(output_mean, output_median, output_mode)
outputTable = signif(outputTable,2)
colnames(outputTable) = NULL
outputTable = as.data.table(outputTable, keep.rownames=T)

outputTable$rn = c('Mean', 'Median', 'Mode')
outputTable$V2 = paste(outputTable$V2, ' (', outputTable$V3, '-', outputTable$V4, ')', sep = '')
outputTable$V5 = paste(outputTable$V5, ' (', outputTable$V6, '-', outputTable$V7, ')', sep = '')
outputTable$V3 = outputTable$V4 = outputTable$V6 = outputTable$V7 = NULL

colnames(outputTable) = c('Summary statistic', 'Observed in B.1.617.2 cases', 'Bootstrap sample mean (95% CI)', 'Differenc (95% CI)')
write.csv(outputTable, 'C:/Users/Sony/Desktop/VOC617/supp_table.csv', row.names = F)
