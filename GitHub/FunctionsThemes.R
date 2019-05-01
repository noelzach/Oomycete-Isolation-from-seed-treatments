# ipak function: install and load multiple R packages.
# check to see if packages are installed. Install them if they are not, then load them into the R session.
# Source: https://gist.github.com/stevenworthington/3178163
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[,"Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

chi.attempt <- function(tissue, year, site){
  attempted_root <- iso[iso$isattempt == TRUE & iso$isloc == as.character(tissue) & iso$Year == year & iso$location ==  as.character(site),]
  
  summ_root_true <- summary(attempted_root$trtname[attempted_root$isrecover == TRUE])
  summ_root_false <- summary(attempted_root$trtname[attempted_root$isrecover == FALSE])
  
  Cmax_success <- summ_root_true[1]
  insuit_success <- summ_root_true[2]
  NTC_success <- summ_root_true[3]
  Cmax_fail <- summ_root_false[1]
  insuit_fail <- summ_root_false[2]
  NTC_fail <- summ_root_false[3]
  
  recov_freq_Cmax_root <- Cmax_success/(Cmax_fail + Cmax_success)
  recov_freq_insuit_root <- insuit_success/(insuit_fail + insuit_success)
  recov_freq_NTC_root <- NTC_success/(NTC_fail + NTC_success)
  
  chi <- matrix(0, nrow = 3, ncol = 2, dimnames = list(c("Cruisermaxx", "IntegoSuite", "NTC"), c("isolation_success", "isolation_failures")))
  chi[,1] <- c(Cmax_success, insuit_success, NTC_success)
  chi[,2] <- c(Cmax_fail, insuit_fail, NTC_fail)
  chi.1 <- addmargins(as.matrix(chi), c(1,2))
  
  
  result <- chisq.test(chi)
  results <- list(chi.1, result)
  return(results)
}

logistic_regression_isolation_sucess <- function(tissue, year, site){
  MI_root_2016 <- iso[iso$isattempt == TRUE & iso$isloc == as.character(tissue) & iso$Year == year & iso$location ==  as.character(site),]
  
  MI_root_2016_logistic <- glm(isrecover ~ trtname, data = MI_root_2016, family = "binomial")
  df_resid <- MI_root_2016_logistic$df.residual
  deviance <- MI_root_2016_logistic$deviance
  overdispersion <- 1-pchisq(deviance,df_resid) # test for overdispersion of the model by testing if the residual deviance is greater than expected based on a Chi-squared distribution
  # it is significant (i.e., < 0.05), meaning there was significant overdispersion
  #hist(residuals(MI_root_2016_logistic)) # check normality of residuals
  #qqnorm(resid(MI_root_2016_logistic)); qqline(resid(MI_root_2016_logistic))
  
  # Test for overall significance of seed treatment on isolation success 
  # Here we are testing if the logistic model with seed treatment as a factor is any better than just fitting the null model with only an intercept.
  LRT <- lmtest::lrtest(MI_root_2016_logistic) 
  
  # Run these if it is significant 
  mod.sum <- summary(MI_root_2016_logistic)
  Tukey <- summary(glht(MI_root_2016_logistic, mcp(trtname="Tukey")))
  return_list <- list(MI_root_2016_logistic, overdispersion, LRT, mod.sum,Tukey)
  return(return_list)
}

my.theme <- theme(axis.text.x = element_text(size = 12, face = "bold", angle=45, hjust=1, family = "serif"),
                  axis.text.y = element_text(size = 15, face = "bold", family = "serif"),
                  axis.title.x = element_text(size = 25, face = "bold", family = "serif"),
                  axis.title.y = element_text(size = 20, face = "bold", family = "serif"),
                  axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
                  axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
                  legend.text = element_text(size = 10, face = "bold", family = "serif"),
                  legend.key = element_blank(),
                  legend.title = element_text(size = 10, face="bold", family = "serif"),
                  legend.position = "right",
                  strip.text.x = element_text(size = 25, face = "bold", family = "serif"),
                  title = element_text(size = 10, family = "serif"))


growth_fun = function(dat) {
  n <- 10
  df <- data.frame(x = dat$hrs, y = dat$od600meanblank)
  spline <- data.frame(spline(df, n=n*10))
  spline$y <- ifelse(spline$y < 0, 0, spline$y)
  ggplot(df, aes(x = x, y = y ) ) +
    geom_line(data=spline, color = "red") + 
    stat_summary(fun.y=mean,geom="point") +
    stat_summary(fun.data = mean_se, geom = "errorbar") +
    stat_summary(fun.y=mean,geom="line") +
    theme_classic() +
    ggtitle(label = paste(dat$is_plot, "Set",dat$set_plot)) + 
    geom_vline(xintercept = 24, linetype = "dashed") +
    geom_vline(xintercept = 48, linetype = "dashed") +
    theme(plot.title = element_text(size = 40, face = "bold"))+
    ylim(c(-0.1, 1.6))
}

extract.zerohour.mean <- function(x) {
  zero.mean <- x %>%
    subset(conc == 0 & chem == "gr") %>%
    select(od600) %>%
    as.numeric()
  return(blank.mean)
}

extract.BLANK.mean <- function(x) {
  blank.mean <- x %>%
    subset(genus == "BLANK") %>%
    select(mean.od600) %>%
    as.numeric()
  return(blank.mean)
}
extract.BLANK.sd <- function(x) {
  sd.BLANK <- x %>%
    subset(genus == "BLANK") %>%
    select(sd.od600) %>%
    as.numeric()
  return(sd.BLANK)
}

Z <- function (sdPC, sdNC, mPC, mNC) {1-((3*(sdPC + sdNC))/abs(mPC - mNC))}

extract.mean.zero <- function(x) {
  mean.0ppm <- x %>%
    subset(conc == 0) %>%
    select(od600meanblank) %>%
    summarise(mean(.$od600meanblank)) %>%
    as.numeric()
  return(mean.0ppm)
}

drm.func <- function(x) {
  try(drm(relgrowth*100 ~ conc, 
          fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), 
          data = x))
}
ED.func <- function(x) {
  estimate <- try(data.frame(ED(x, type = "absolute", respLev = 50, display = FALSE)))
  EC50 <- try(estimate[[1]])
  return(EC50)
}
drm.func.convergence <- function(x){
  try(drm(relgrowth*100 ~ conc, 
      fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), 
      data = x, start = c(0.5, 0, 100, 0.01)))
}
n_fun <- function(x){return(data.frame(y = mean(x)*0 - 0.03, label = paste0(length(x))))}

