# code for cluster-based permutation analysis
library(dplyr)
library(arm)
# prep data and define timeline
ddata <- ddata[ddata$condition2=="critical" & ddata$time>=-1000 & ddata$time<=1000,]
data <- group_by(ddata, language.group, subject, condition4, time) %>% summarise(proximallook.mean = mean(proximallook))
data$condition <- data$condition4
data$proximallook<-ifelse(data$proximallook.mean==0,logit(0.025),ifelse(data$proximallook.mean==1,logit(0.975),logit(data$proximallook.mean)))
StartTime = -1000
EndTime = 1000
StepSize = 100
# function to calculate *observed* values
calculate.statistics = function(data,StartTime,EndTime,StepSize){
  statsCol<- rep(NA,length(seq(StartTime,EndTime,StepSize)))
  # create a table to store each regression result
  stats1 = data.frame(time = statsCol, Estimate= statsCol,SE = statsCol, Stat = statsCol, pval = statsCol)
  # j is the order, i is the actual timestamp, stepSize maps i to j
  j = 1
  for (i in seq(StartTime,EndTime,StepSize)){
    # run regression
    dat.glm = lm(proximallook~condition, data = data, subset = time ==(i))
    # save regression results
    stats1[j,1:4] = c(i, coef(summary(dat.glm))[2,1:3])
    # set t-value threshold for clusters to 1.6 (as in Wittenberg, Khan & Snedeker, 2017)
    stats1[j,5] = ifelse(abs(stats1[j,4]) > 1.6, 0.04, 0.5)
    j = j+1
  }
  return(list(stats1))
}
# function to calculate *permuted* values
resample = function(data,StartTime,EndTime,StepSize){
  for (i in unique(data$subject)){
    # subset to each subject
    data_S <- data[data$subject == i,]
    # randomly sample (permute) condition labels
    data_S$condition <- sample(data_S$condition, length(data_S$condition))
    data[data$subject == i,] <- data_S
  }
  data$condition = as.factor(data$condition)
  # run regression with permuted data
  stats.resample = calculate.statistics(data,StartTime,EndTime,StepSize)
  return(stats.resample)
}
# function to identify clusters
find.clusters = function(pval.vector, tval.vector, latencies, alpha=.05, cluster.size=2,signed=FALSE) {
  binary.stat <- as.numeric(pval.vector < alpha)
  tval.stat <- rep(0,length(binary.stat))
  tval.stat[2:length(binary.stat)] = ((tval.vector[2:length(binary.stat)] * tval.vector[1:(length(binary.stat)-1)]) <0)
  cluster.start <- c()
  cluster.end <- c()
  cluster.stat <- c()
  in.cluster <- 0
  found.clusters <- 0
  new.end <- 0
  new.start <- 0
  end.cluster <- 0
  for (n in 1:length(binary.stat)) {
    if (signed == FALSE){
      if (in.cluster == 0 && binary.stat[n] == 1 ) {
        new.start = n
        in.cluster = 1
      }
    }else{
      if (in.cluster == 0 && binary.stat[n] == 1 && tval.vector[n] >= 0) {
        new.start = n
        in.cluster = 1
      }
    }
    if (in.cluster == 1 && binary.stat[n] == 0) {
      new.end = n
      end.cluster = 1
      in.cluster = 0
    }
    if (in.cluster == 1 && tval.stat[n] == 1) {
      new.end = n
      end.cluster = 1
      in.cluster = 0
    }
    if (in.cluster == 1 && binary.stat[n] == 1 && n == length(binary.stat)) {
      new.end = n
      end.cluster = 1
      in.cluster = 1
    }
    if (end.cluster) {
      if ((new.end - new.start) >= cluster.size) {
        found.clusters <- found.clusters + 1
        cluster.start<- c(cluster.start, latencies[new.start])
        cluster.end <- c(cluster.end, latencies[new.end])
        cluster.stat <- c(cluster.stat, sum(abs(tval.vector[new.start:(new.end)])))
      }
      end.cluster = 0
    }
  }
  cluster.out <- data.frame(start = cluster.start, end = cluster.end, cluster.stat = cluster.stat)
  return(cluster.out)
}
# function to obtain p-values
pval = function(sample.cluster,resample.large.cluster){
  sample.cluster$pval = 1
  for (i in 1:length(sample.cluster$cluster.stat)){
    sample.cluster$pval[i] = 1 - sum(as.numeric(sample.cluster$cluster.stat[i]>resample.large.cluster)/Resamples_N)
  }
  return(sample.cluster)
}
# apply the above functions for each language.group
x <- unique(data$language.group)
for (a in 1:length(x)) {
  group.data <- subset(data, language.group == x[a])
  # calculate *observed* values and clusters
  stats.sample = calculate.statistics(group.data,StartTime,EndTime,StepSize)
  sample.cluster1 = find.clusters(stats.sample[[1]]$pval,stats.sample[[1]]$Stat,stats.sample[[1]]$time,cluster.size=2)
  # calculate *permuted* values and clusters
  Resamples_N = 1000
  resample1.large.cluster = rep(NA,Resamples_N)
  for (i in 1:Resamples_N){
    stats.resample = resample(data,StartTime,EndTime,StepSize)
    resample1.max = ifelse(is.numeric(find.clusters(stats.resample[[1]]$pval,stats.resample[[1]]$Stat,stats.resample[[1]]$time,cluster.size=2)$cluster.stat),max(find.clusters(stats.resample[[1]]$pval,stats.resample[[1]]$Stat,stats.resample[[1]]$time,cluster.size=2)$cluster.stat),0)
    # find the largest cluster (are observed clusters > 95% of largest clusters?)
    resample1.large.cluster[i] = ifelse(resample1.max>0,resample1.max,0)
  }
  # calculate p-values (if there were any significant observed clusters)
  if (length(sample.cluster1$cluster.stat)>=1) {
    sample.cluster1 = pval(sample.cluster1,resample1.large.cluster)
    sample.cluster1$pval <- ifelse(sample.cluster1$pval==0, "< 0.001", sample.cluster1$pval)
    filename <- paste((x[a]), ".proximity.clusters.csv", sep="")
    write.csv(sample.cluster1, file=filename, row.names=F)
  }
}