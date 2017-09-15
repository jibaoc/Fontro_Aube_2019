##################################################
# The goal of this script is :
#   - to fit the proportion of a feature in each up-regulated exons to a zero-inflated beta model
#   - to fit the proportion of a feature in each down-regulated exons to a zero-inflated beta model
#   - compaire those to models
#   - if they are different, then the up and down regulated exons show different distribution of proportions for this feature



###################################################
# Getting arguments
args = commandArgs(trailingOnly = TRUE)
if(length(args) > 1){
  res = as.vector(args)
}else{
  stop("Error")
}
########################################################
# Preparation of the vectors for comparison
feature_up = as.numeric(strsplit(res[1],",")[[1]])

feature_down = as.numeric(strsplit(res[2],",")[[1]])


feature_all = c(feature_up, feature_down)
feature_group = c(rep("up", times = length(feature_up)), rep("down", times = length(feature_down)))



dataf = data.frame(cbind(feature_all, feature_group))
dataf$feature_all = as.numeric(as.vector(dataf$feature_all))

##############################################################
library(gamlss)
fit1 = gamlss(feature_all ~ 1, family=BEINF0, data=dataf)
fit2 = gamlss( feature_all ~ feature_group, sigma.formula=~feature_group, nu.formula = ~feature_group, family=BEINF0, data=dataf)

cat(LR.test(fit1, fit2))


nu_up = as.numeric(fitted(fit2,"nu")[1])
mu_up = as.numeric(fitted(fit2,"mu")[1])
sigma_up = as.numeric(fitted(fit2,"sigma")[1])


nu_down = as.numeric(fitted(fit2,"nu")[length(fitted(fit2,"nu"))])
mu_down = as.numeric(fitted(fit2,"mu")[length(fitted(fit2,"mu"))])
sigma_down = as.numeric(fitted(fit2,"sigma")[length(fitted(fit2,"sigma"))])




cat("sigma_up=")
cat(sigma_up)
cat("\nmu_up=")
cat(mu_up)
cat("\nnu_up=")
cat(nu_up)
cat("\nsigma_down=")
cat(sigma_down)
cat("\nmu_down=")
cat(mu_down)
cat("\nnu_down=")
cat(nu_down)