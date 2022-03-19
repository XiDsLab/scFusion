## README ####
### V 2.2.1
##############

options(warn=-1)
## 设分布都是ZINB分布，先用极大似然估计参数
CalcZINB = function(par)
{
  sum = 0
  for (i in 1:dim(data_filtered)[1])
  {
    cc = as.numeric(stringr::str_split(stringr::str_sub(data_filtered$V4[i], 2, -2), ', ')[[1]])
    expr1 = as.numeric(stringr::str_split(stringr::str_sub(data_filtered$V15[i], 2, -2), ', ')[[1]])
    expr2 = as.numeric(stringr::str_split(stringr::str_sub(data_filtered$V16[i], 2, -2), ', ')[[1]])
    if (sum(is.na(expr1))>0 | sum(is.na(expr2))>0 | sum(is.na(cc))>0){
      next()
    }
    aveexpr1 = mean(expr1)
    aveexpr2 = mean(expr2)
    gc = (data_filtered$V13[i] + data_filtered$V14[i]) / 2
    gcpre = log(gcprelist[round(gc * 1000) + 1])
    lambda = exp(par[5])
    linpart = par[6] + par[7] * gcpre + par[8] * aveexpr1 + par[9] * aveexpr2
    zeroprob = min(exp(linpart) / (1 + exp(linpart)), 0.9999999999999)
    for (j in 1:length(cc)) {
      mu = par[1] + par[2] * gcpre + exp(par[3]) * expr1[j] + exp(par[4]) * expr2[j]
      mu = exp(mu)
      sum = sum + dnbinom(as.numeric(cc[j]),
                          1 / lambda,
                          prob = 1 / (1 + mu * lambda),
                          log = TRUE) + log(1 - zeroprob)
      sum = sum - log(1 - dnbinom(0, 1 / lambda, prob = 1 / (1 + mu * lambda))) # 因为是truncated NB
    }
    sum = sum + (numcell - length(cc)) * log(zeroprob)
    sum = sum - log(1 - zeroprob ^ numcell - numcell * (1 - zeroprob) *
                      zeroprob ^ (numcell - 1)) #条件在至少两个细胞支持
    if (is.nan(sum)){
      break
    }
  }
  - sum
}


CalcPValueZINB_Fusion_2.0.0_3 = function(par)
{
  Pv = rep(0, dim(data)[1])
  for (i in 1:dim(data)[1]){
    cc = as.numeric(stringr::str_split(stringr::str_sub(data$V4[i], 2, -2), ', ')[[1]])
    if (sum(cc) <= 3 | data$V3[i] <= 1){
      Pv[i] = 1
      next()
    }
    expr1 = as.numeric(stringr::str_split(stringr::str_sub(data$V15[i], 2, -2), ', ')[[1]])
    expr2 = as.numeric(stringr::str_split(stringr::str_sub(data$V16[i], 2, -2), ', ')[[1]])
    if (sum(is.na(expr1))>0 | sum(is.na(expr2))>0 | sum(is.na(cc))>0){
      Pv[i] = 1
      next()
    }
    aveexpr1 = mean(expr1) / 2
    aveexpr2 = mean(expr2) / 2
    gc = (data$V13[i] + data$V14[i]) / 2
    gcpre = log(gcprelist[round(gc * 1000) + 1])
    numcellsup = length(cc)
    linpart = par[6] + par[7] * gcpre + par[8] * aveexpr1 + par[9] * aveexpr2
    zeroprob = min(exp(linpart) / (1 + exp(linpart)), 0.9995)
    betamu = 1 - zeroprob
    betavar = betamu * (1 - betamu) / numcell
    myalpha = (betamu * (1 - betamu) / betavar - 1) * betamu
    mybeta = (1 / betamu - 1) * myalpha
    for (j in 1:length(cc)) {
      cc[j] = min(15, cc[j])
    }
    totalread = sum(cc)
    mu = par[1] + par[2] * gcpre + exp(par[3]) * 2 * aveexpr1 + exp(par[4]) * 2 * aveexpr2
    truemu = exp(mu)
    mu = max(truemu, 0.1)
    lambda = exp(par[5])
    truevar = truemu + truemu ^ 2 * lambda
    lambda = max(0.001, (truevar - mu) / (mu ^ 2))
    simuset = rnbinom(10000, 1 / lambda, prob = 1 / (1 + mu * lambda))
    simuset = simuset[simuset > 0]
    whilecount = 0
    while (length(simuset) < 2000){
      newsimuset = rnbinom(10000, 1 / lambda, prob = 1 / (1 + mu * lambda))
      newsimuset = newsimuset[newsimuset>0]
      simuset = c(simuset, newsimuset)
      whilecount = whilecount + 1
      if (whilecount > 100){
        simuset = c(simuset, 1)
      }
    }
    simuset[simuset >= 15] = 15
    totalreadsample = c()
    for (j in 1:100){
      nonzeronum = rbinom(1, numcell, 1-zeroprob)
      usedreadnum = sample(simuset, nonzeronum, replace = T)
      totalreadsample = c(totalreadsample, sum(usedreadnum))
    }
    totalreadmean = mean(totalreadsample)
    totalreadsd = sd(totalreadsample)
    thispvalue = pnorm(totalread, totalreadmean, totalreadsd, lower.tail = F)
    Pv[i] = thispvalue
  }
  #Pv = p.adjust(Pv, method = 'fdr')
  Pv
}

library(splines)
library(stringr)
set.seed(1122)
Args = commandArgs()
numcell = as.integer(Args[6])
data = read.table(Args[7],sep='\t', stringsAsFactors = F)
data_filtered = read.table(Args[8],sep='\t', stringsAsFactors = F)
lastpar = rep(-100, 9)
if (length(Args) >= 10){
  if (file.exists(Args[10])){
    load(Args[10])
    lastpar = optres$par
  }
}

### filter bad gc
neiflag = c()
for (i in 1:dim(data)[1]) {
  if (max(data$V14[i], data$V13[i]) <= 0.05||min(data$V14[i], data$V13[i]) >= 0.95) {
    neiflag[i] = 1
  } else{
    neiflag[i] = 0
  }
}
data = data.frame(data, neiflag)
data = data[data$neiflag == 0,]
row.names(data) = 1:dim(data)[1]


# fit a spline for gc
gc = round((data$V14 + data$V13) / 2, digits = 3)
gclist = sort(unique(gc))
gclist = round(gclist, digits = 3)
gclist = unique(gclist)
chimexpr = gclist
for (i in 1:length(gclist)) {
  usedata = data[gc == gclist[i],]
  exprtotal = 0
  count = 0
  for (j in 1:dim(usedata)[1]) {
    cc = as.numeric(stringr::str_split(stringr::str_sub(usedata$V4[j], 2, -2), ', ')[[1]])
    if (length(cc) >= 1) {
      exprtotal = exprtotal + sum(cc)
      count = count + length(cc)
    }
  }
  chimexpr[i] = exprtotal / count
}
lmres = lm(chimexpr ~ bs(gclist, df = 5, intercept = T))
basepre = predict(lmres, data.frame(gclist = gclist))
diff = abs(chimexpr - basepre)
thred = mean(diff) + 1.5 * sqrt(var(diff))
for (i in 1:length(gclist)) {
  pre = basepre[i]
  if (chimexpr[i] - pre > thred) {
    chimexpr[i] = pre
  }
}
lmres = lm(chimexpr ~ bs(gclist, df = 5, intercept = T))
gcprelist = predict(lmres, data.frame(gclist = seq(0, 1, length.out = 1001)))
gcprelist[gcprelist <= 0] = 0.01

print('Start Reading ChiDist File!')
neiflag = rep(0, dim(data)[1])
totalcount = rep(0, dim(data)[1])
datacc = rep('', dim(data)[1])
dataexpmin = rep('', dim(data)[1])
dataexpmax = rep('', dim(data)[1])
datagcmin = rep(0, dim(data)[1])
datagcmax = rep(0, dim(data)[1])
for (i in 1:dim(data)[1]) {
  cc = as.numeric(stringr::str_split(stringr::str_sub(data$V4[i], 2, -2), ', ')[[1]])
  totalcount[i] = sum(cc)
  smallindex = which(cc <= 2)
  bigindex = which(cc > 2)
  if (max(cc) < 2 | length(cc) < 2) {
    neiflag[i] = 1
  }
  count1 = sum(cc<=2)
  usesmallcellnum = floor(data$V5[i] * 1.1 + 2)
  usesmallcell = c()
  if (count1 <= usesmallcellnum){
    usesmallcell = smallindex
  }else{
    usesmallcell = sample(smallindex, usesmallcellnum)
  }
  cc = cc[union(bigindex, usesmallcell)]
  if (neiflag[i] == 1) {
    next()
  }
  expmin = c()
  expmax = c()
  gene1expr = as.numeric(stringr::str_split(stringr::str_sub(data$V15[i], 2, -2), ', ')[[1]])
  gene2expr = as.numeric(stringr::str_split(stringr::str_sub(data$V16[i], 2, -2), ', ')[[1]])
  if (min(min(gene1expr), min(gene2expr)) < 0.0001) {
    neiflag[i] = 1
    next()
  }
  datagcmin[i] = min(data$V14[i], data$V13[i])
  datagcmax[i] = max(data$V14[i], data$V13[i])
  for (j in 1:length(gene1expr)) {
    expmin = c(expmin, log(1 + min(
      as.numeric(gene1expr[j]), as.numeric(gene2expr[j])
    ) / 100))
    expmax = c(expmax, log(1 + max(
      as.numeric(gene1expr[j]), as.numeric(gene2expr[j])
    ) / 100))
  }
  expmin = expmin[union(bigindex, usesmallcell)]
  expmax = expmax[union(bigindex, usesmallcell)]
  datacc[i] = paste0('[', str_c(as.character(cc), collapse = ', '), ']')
  dataexpmin[i] = paste0('[', str_c(as.character(expmin), collapse = ', '), ']')
  dataexpmax[i] = paste0('[', str_c(as.character(expmax), collapse = ', '), ']')
}
data$V4 = datacc
data$V15 = dataexpmin
data$V16 = dataexpmax
data$V13 = datagcmin
data$V14 = datagcmax
data$neiflag = neiflag
data$totalreadcount = totalcount
data = data[data$neiflag == 0,]
row.names(data) = 1:dim(data)[1]

if (norm(lastpar - rep(-100, 9), '2') <= 1){
  neiflag = rep(0, dim(data_filtered)[1])
  totalsup = rep(0, dim(data_filtered)[1])
  lastpos1 = c(-1)
  lastpos2 = c(-1)
  cellsupqs = quantile(data_filtered$V3)
  highthres = cellsupqs[3] + cellsupqs[4] - cellsupqs[2]
  for (i in 1:dim(data_filtered)[1]) {
    cc = as.numeric(stringr::str_split(stringr::str_sub(data_filtered$V4[i], 2, -2), ',')[[1]])
    m = mean(cc)
    a = runif(1)
    if (a > m ^ 0.1 * 70000 / dim(data_filtered)[1]) {
      neiflag[i] = 1
      next()
    }
    if (max(cc) < 2 && dim(data_filtered)[1]>150) {
      neiflag[i] = 0
    }else{
      neiflag[i] = 0
    }
    if (max(cc) < 2 && length(cc) <= 5){
      neiflag[i] = 1
    }
    if (data_filtered$V3[i] > highthres && dim(data_filtered)[1]<150){
      neiflag[i] = 1
    }
    for (k in 1:length(lastpos1)) {
      if (data_filtered$V9[i] == lastpos1[k] &&
          data_filtered$V10[i] == lastpos2[k]||data_filtered$V9[i] == lastpos2[k] &&
          data_filtered$V10[i] == lastpos1[k]) {
        neiflag[i] = 1
        break()
      }
    }
    if (neiflag[i] == 0) {
      lastpos1 = c(lastpos1, data_filtered$V9[i])
      lastpos2 = c(lastpos2, data_filtered$V10[i])
      if (length(lastpos1) > 100) {
        lastpos1 = lastpos1[2:101]
        lastpos2 = lastpos2[2:101]
      }
    }
    if (neiflag[i] == 1) {
      next()
    }
    totalsup[i] = sum(cc)
    data_filtered$V13[i] = min(data_filtered$V14[i], data_filtered$V13[i])
    data_filtered$V14[i] = max(data_filtered$V14[i], data_filtered$V13[i])
    expmin = c()
    expmax = c()
    gene1expr = as.numeric(stringr::str_split(stringr::str_sub(data_filtered$V15[i], 2, -2), ', ')[[1]])
    gene2expr = as.numeric(stringr::str_split(stringr::str_sub(data_filtered$V16[i], 2, -2), ', ')[[1]])
    for (j in 1:length(gene1expr)) {
      expmin = c(expmin, log(1 + min(
        as.numeric(gene1expr[j]), as.numeric(gene2expr[j])
      ) / 100))
      expmax = c(expmax, log(1 + max(
        as.numeric(gene1expr[j]), as.numeric(gene2expr[j])
      ) / 100))
    }
    data_filtered$V15[i] = paste0('[', str_c(as.character(expmin), collapse = ', '), ']')
    data_filtered$V16[i] = paste0('[', str_c(as.character(expmax), collapse = ', '), ']')
    if (data_filtered$V7[i] == data_filtered$V8[i] &
        abs(data_filtered$V10[i] - data_filtered$V9[i]) < 200000) {
      neiflag[i] = 1
    } else{
      neiflag[i] = 0
    }
    if (min(min(gene1expr), min(gene2expr)) < 0.0001) {
      neiflag[i] = 1
    }
  }
  data_filtered = data.frame(data_filtered, neiflag)
  data_filtered = data_filtered[data_filtered$neiflag == 0,]
  row.names(data_filtered) = 1:dim(data_filtered)[1]
  
  neiflag = rep(0, dim(data_filtered)[1])
  for (i in 1:dim(data_filtered)[1]) {
    cc = as.numeric(stringr::str_split(stringr::str_sub(data_filtered$V4[i], 2, -2), ',')[[1]])
    m = mean(cc)
    a = runif(1)
    if (a > m ^ 0.1 * 700 / dim(data_filtered)[1]) {
      neiflag[i] = 1
    } else{
      neiflag[i] = 0
    }
  }
  data_filtered$neiflag = neiflag
  data_filtered = data_filtered[data_filtered$neiflag == 0,]
  row.names(data_filtered) = 1:dim(data_filtered)[1]
}



if (norm(lastpar - rep(-100, 9), '2') < 1){
  print('Start Estimating Parameters!')
  optres = optim(rnorm(9, mean=-0.4, sd=0.1), CalcZINB, control = list(maxit = 10000))
  if (length(Args) >= 10){
    save(optres, file = Args[10])
  }
}else{
  print('Skip Estimating Parameters!')
}


print('Start Calculating P-values!')
PvalueFusion23 = CalcPValueZINB_Fusion_2.0.0_3(optres$par)
data = cbind(data, PvalueFusion23)
TrueFusion = rep(1, dim(data)[1])

print('Preparing for Output!')
goodindex = which(TrueFusion <= 2)
totalnumfusion = length(goodindex)
allresult = data.frame(
  'FusionName' = rep(NA, totalnumfusion),
  CellSupport = rep(NA, totalnumfusion),
  JunctionReadCount = rep(NA, totalnumfusion),
  SpanningFragCount = rep(NA, totalnumfusion),
  Position1 = rep(NA, totalnumfusion),
  Position2 = rep(NA, totalnumfusion),
  FakeProb = rep(NA, totalnumfusion),
  Pv23 = rep(NA, totalnumfusion),
  strand1 = rep(NA, totalnumfusion),
  strand2 = rep(NA, totalnumfusion)
)

fusioncandidateread = data.frame(
  JunctionRead = rep(NA, totalnumfusion),
  Brkpnt = rep(NA, totalnumfusion),
  Prob = rep(NA, totalnumfusion)
)
fusioncandidate = rep('', totalnumfusion)
for (i in 1:totalnumfusion) {
  fusioncandidate[i] = paste(data$V1[i], data$V2[i], sep = '--')
}
allresult$FusionName = fusioncandidate
allresult$JunctionReadCount = data$totalreadcount
fusioncandidateread$JunctionRead = data$V19
fusioncandidateread$Brkpnt = data$V20
fusioncandidateread$Prob = data$V25
allresult$FakeProb = data$V25
allresult$Pv23 = data$PvalueFusion23
allresult$CellSupport = data$V3
allresult$strand1 = data$V21
allresult$strand2 = data$V22
spanningfragcount = rep(0, totalnumfusion)
pos1 = rep('', totalnumfusion)
pos2 = rep('', totalnumfusion)
for (i in 1:totalnumfusion) {
  if (as.numeric(data$V5[i]) == 0) {
    spanningfragcount[i] = 0
  } else{
    spanningfragcount[i] = sum(as.numeric(stringr::str_split(
      stringr::str_sub(data$V6[i], 2, -2), ', '
    )[[1]]))
  }
  pos1[i] = paste(data$V7[i], data$V9[i], sep = ':')
  pos2[i] = paste(data$V8[i], data$V10[i], sep = ':')
}
allresult$SpanningFragCount = spanningfragcount
allresult$Position1 = pos1
allresult$Position2 = pos2
write.table(allresult, Args[9], sep = '\t', quote = F, row.names = F)
