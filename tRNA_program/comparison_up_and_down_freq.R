###################################################
# Getting arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args) > 1){
  res = as.vector(args)
}else{
  stop("Error")
}

feature_up = as.numeric(strsplit(res[1],",")[[1]])

feature_down = as.numeric(strsplit(res[2],",")[[1]])

###################################################

feature_all = c(feature_up, feature_down)
feature_group = c(rep("up", times = length(feature_up)), rep("down", times = length(feature_down)))

library(gamlss)

dataf = data.frame(cbind(feature_all, feature_group))
dataf$feature_all = as.numeric(as.vector(dataf$feature_all))

fit1 = gamlss(feature_all ~ 1, family=BEINF0, data=dataf)
fit2 = gamlss( feature_all ~ feature_group, sigma.formula=~feature_group, nu.formula = ~feature_group, family=BEINF0, data=dataf)

cat(LR.test(fit1, fit2))
