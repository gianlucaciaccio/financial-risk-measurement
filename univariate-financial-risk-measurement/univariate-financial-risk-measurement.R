library(readxl)
library(zoo)
library(quarks)
library(rugarch)
library(moments)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
source("functions.R")


#### Load the data and calculate returns ####
sp500close <- read_excel("data.xlsx",col_types = c("date","numeric"))

sp500close$date <- as.Date(sp500close$date)

LogReturn <- diff(log(sp500close$close))

sp500date <- sp500close$date[-1]


#### Plot: histogram log returns and Normal density ####
data.frame(LogReturn) %>%
  ggplot(mapping = aes(LogReturn)) +
  geom_histogram(aes(y = ..density..), bins = 50,
                 fill = "cornflowerblue", color="gray15") +
  stat_function(fun = dnorm, color="red",
                args = list(mean = mean(LogReturn),
                            sd = sd(LogReturn))) +
  theme_gianluca() +
  labs(y="",title = "Distribution of Log Returns vs. Normal Density")

#### Plot: Q-Q Plot log returns and Normal quantiles ####
ggplot(mapping = aes(sample=LogReturn)) +
  geom_qq(color="cornflowerblue", shape=1) + 
  geom_qq_line(color="red") +
  theme_gianluca() +
  labs(x="Theorical", y="Sample",
       title = "Normal Q-Q Plot of Log Returns")



#### Calculate moments of log return and perform test for Normality ####
sp500close %>%
  mutate(return = log(close/lag(close))) %>%
  filter(!is.na(return)) %>%
  summarize(mean = mean(return), sd = sd(return), 
            skewness = skewness(return), 
            kurtosis = kurtosis(return),
            EK = kurtosis(return)-3)

jarque.test(LogReturn)



#### Check Student's distribution of log returns ####
data_qqt <- data.frame(
  LogReturn=qqt(LogReturn,plot.it = FALSE)$y,
  QQt_1=qqt(LogReturn, df = 1, plot.it = FALSE)$x,
  QQt_2=qqt(LogReturn, df = 2, plot.it = FALSE)$x,
  QQt_3=qqt(LogReturn, df = 3, plot.it = FALSE)$x,
  QQt_4=qqt(LogReturn, df = 4, plot.it = FALSE)$x,
  QQt_5=qqt(LogReturn, df = 5, plot.it = FALSE)$x,
  QQt_6=qqt(LogReturn, df = 6, plot.it = FALSE)$x,
  QQt_7=qqt(LogReturn, df = 7, plot.it = FALSE)$x,
  QQt_8=qqt(LogReturn, df = 8, plot.it = FALSE)$x,
  QQt_9=qqt(LogReturn, df = 9, plot.it = FALSE)$x,
  QQt_10=qqt(LogReturn, df = 10, plot.it = FALSE)$x) %>%
  pivot_longer(cols = -c("LogReturn"), names_to = c("QQt","df"), 
               names_sep = "_", values_to = "Student_QT") %>%
  select(LogReturn,df,Student_QT) %>%
  mutate(df = as.factor(paste(df,"df",sep = " ")), 
         df = fct_relevel(df,"10 df",after = Inf))


data_qqtline <- data_qqt %>%
  group_by(df) %>%
  summarize(qt75 = qt(c(0.75),df = as.numeric(df)),
            qt25 = qt(c(0.25),df = as.numeric(df)),
            quant75 = quantile(LogReturn,probs=0.75),
            quant25 = quantile(LogReturn,probs=0.25),
            slope = (quant75-quant25)/(qt75-qt25),
            intercept = quant25-slope*qt25) %>%
  distinct() %>%
  select(df,slope,intercept) %>%
  mutate()

#### Plot: Q-Q Plot log returns and t quantiles by different d.f. ####
ggplot(data_qqt, aes(x=Student_QT,y=LogReturn)) +
  geom_point(shape=1, color = "cornflowerblue") +
  geom_abline(data=data_qqtline,
              aes(slope=slope,intercept=intercept),
              colour = "red", size = 0.55) +
  facet_wrap(~df, scales = "free") +
  theme_gianluca() +
  theme(panel.grid = element_blank()) +
  labs(x="Theorical", y="Sample",
       title = "Student's Q-Q Plots of Log Returns by different d.f.")


#### Plot: Time series of SP500 closing price and log returns ####
sp500close %>%
  mutate(return = log(close/lag(close))) %>%
  pivot_longer(cols = -c("date"), names_to = "Series",
               values_to = "Value") %>%
  mutate(Series = ifelse(Series=="close","Closing Price","Log Return")) %>%
  ggplot(aes(x = date, y = Value)) +
  geom_line(color = "cornflowerblue", size = 0.55) +
  facet_wrap(~Series, ncol = 1, scales = "free") +
  theme_gianluca() +
  labs(x = "", y = "", title = "SP500 time series 2001-2010: Closing price and log return")


#### ACF of log returns and squared log returns ####
data_autocor <- data.frame(
  Lag = acf(LogReturn,lag.max = 100,plot = FALSE)$lag[-1],
  Log_Return = acf(LogReturn,lag.max = 100,plot = FALSE)$acf[-1],
  Squared_Log_Return = acf(LogReturn^2,lag.max = 100,plot = FALSE)$acf[-1]) %>%
  pivot_longer(cols = -c("Lag"), names_to = "Series", values_to = "ACF") %>%
  mutate(Series = chartr("_", " ", Series))

#### Plot: ACF of log returns and squared log returns ####
ggplot(data_autocor, aes(x=Lag,y=ACF)) +
  geom_line(color = "cornflowerblue", size = 0.55) +
  facet_wrap(~Series, ncol = 1, scales = "free") +
  theme_gianluca() +
  labs(title = "Autocorrelations")


#### Calculate VaR and ES series using rolling windows ####
data_window <- data.frame(
  Day = sp500date[1000:2513],LogReturn = LogReturn[1000:2513],
  HS_250_VaR = rollapply(LogReturn,width=250, FUN = function(x) 
    -hs(x, p = 0.99, method = "plain")$VaR)[-c(1:750)],
  HS_250_ES = rollapply(LogReturn,width=250, FUN = function(x) 
    -hs(x, p = 0.99, method = "plain")$ES)[-c(1:750)],
  HS_1000_VaR = rollapply(LogReturn,width=1000, FUN = function(x) 
    -hs(x, p = 0.99, method = "plain")$VaR),
  HS_1000_ES = rollapply(LogReturn,width=1000, FUN = function(x) 
    -hs(x, p = 0.99, method = "plain")$ES),
  WHS_250_VaR = rollapply(LogReturn,width=250, FUN = function(x) 
    -hs(x, p = 0.99, method = "age", lambda = 0.95)$VaR)[-c(1:750)],
  WHS_250_ES = rollapply(LogReturn,width=250, FUN = function(x) 
    -hs(x, p = 0.99, method = "age", lambda = 0.95)$ES)[-c(1:750)],
  WHS_1000_VaR = rollapply(LogReturn,width=1000, FUN = function(x) 
    -hs(x, p = 0.99, method = "age", lambda = 0.95)$VaR),
  WHS_1000_ES = rollapply(LogReturn,width=1000, FUN = function(x) 
    -hs(x, p = 0.99, method = "age", lambda = 0.95)$ES),
  Norm_250_VaR = rollapply(LogReturn,width=250, FUN = function(x) 
    uVaR_ES(y=x, type = "norm", p = 0.01)$VaR)[-c(1:750)],
  Norm_250_ES = rollapply(LogReturn,width=250, FUN = function(x) 
    uVaR_ES(y=x, type = "norm", p = 0.01)$ES)[-c(1:750)],
  Norm_1000_VaR = rollapply(LogReturn,width=1000, FUN = function(x) 
    uVaR_ES(y=x, type = "norm", p = 0.01)$VaR),
  Norm_1000_ES = rollapply(LogReturn,width=1000, FUN = function(x) 
    uVaR_ES(y=x, type = "norm", p = 0.01)$ES),
  t_250_VaR = rollapply(LogReturn,width=250, FUN = function(x) 
    uVaR_ES(y=x, type = "t", p = 0.01, v = 3)$VaR)[-c(1:750)],
  t_250_ES = rollapply(LogReturn,width=250, FUN = function(x) 
    uVaR_ES(y=x, type = "t", p = 0.01, v = 3)$ES)[-c(1:750)],
  t_1000_VaR = rollapply(LogReturn,width=1000, FUN = function(x) 
    uVaR_ES(y=x, type = "t", p = 0.01, v = 3)$VaR),
  t_1000_ES = rollapply(LogReturn,width=1000, FUN = function(x) 
    uVaR_ES(y=x, type = "t", p = 0.01, v = 3)$ES),
  CF_250_VaR = rollapply(LogReturn,width=250, FUN = function(x) 
    uVaR_ES(y=x, type = "cf", p = 0.01)$VaR)[-c(1:750)],
  CF_250_ES = rollapply(LogReturn,width=250, FUN = function(x) 
    uVaR_ES(y=x, type = "cf", p = 0.01)$ES)[-c(1:750)],
  CF_1000_VaR = rollapply(LogReturn,width=1000, FUN = function(x) 
    uVaR_ES(y=x, type = "cf", p = 0.01)$VaR),
  CF_1000_ES = rollapply(LogReturn,width=1000, FUN = function(x) 
    uVaR_ES(y=x, type = "cf", p = 0.01)$ES)) %>% 
  pivot_longer(cols = -c("Day","LogReturn"),
               names_to = c("Method","Window","Measure"),
               names_sep = "_", values_to = "Value") %>%
  mutate_if(is.character, as.factor) %>%
  mutate(Method = fct_relevel(Method,c("HS","WHS","Norm","t","CF")),
         Window = fct_relevel(Window,c("250","1000")),
         Window = fct_recode(Window,"250 days"="250","1000 days"="1000"),
         Measure = fct_relevel(Measure,c("VaR","ES")))

#### Plot: VaR violation HS and WHS ####
data_window %>%
  mutate(Violation = ifelse(LogReturn<Value,"VaR violation",
                            "Daily returns SP500")) %>%
  ggplot() +
  geom_point(aes(Day,LogReturn,fill=Violation), shape =21, alpha=0.3) +
  scale_fill_manual(values = c("gray50","red"))+
  geom_line(aes(Day,Value, color = Measure), size = 0.6) +
  facet_grid(Window~Method) + 
  coord_cartesian(y=c(-0.11, -0.002)) +
  scale_color_viridis_d(option = "plasma")+
  labs(y="", fill="",
       title = "VaR violations (p = 0.01) on rolling windows") +
  theme_gianluca()





#### Calculate VaR by different p (from 0.01 to 0.1) ####
data_threshold <- vector("list",100)

for (p in seq(0.001,0.1,0.001)) {
  data_threshold[[p*1000]] <- data.frame(
    p = p,
    HS_VaR = -hs(LogReturn, p = 1-p, method = "plain")$VaR,
    HS_ES = -hs(LogReturn,p=1-p,method = "plain")$ES,
    WHS_VaR = -hs(LogReturn,p = 1-p, method = "age",lambda = 0.95)$VaR,
    WHS_ES = -hs(LogReturn,p=1-p, method = "age",lambda = 0.95)$ES,
    Norm_VaR = uVaR_ES(y = LogReturn, type = "norm", p = p)$VaR,
    Norm_ES = uVaR_ES(y = LogReturn, type = "norm", p = p)$ES,
    t_VaR = uVaR_ES(y = LogReturn, type = "t", p = p, v = 3)$VaR,
    t_ES = uVaR_ES(y = LogReturn, type = "t", p = p, v = 3)$ES,
    CF_VaR = uVaR_ES(y = LogReturn, type = "cf", p = p)$VaR,
    CF_ES = uVaR_ES(y = LogReturn, type = "cf", p = p)$ES)
}

data_threshold <- do.call(rbind,data_threshold) %>%
  pivot_longer(cols = -c("p"), names_to = c("Method","Measure"),
               names_sep = "_", values_to = "Value") %>%
  mutate_if(is.character, as.factor) %>%
  mutate(Method = fct_relevel(Method,c("HS","WHS","Norm","t","CF")),
         Measure = fct_relevel(Measure,c("VaR","ES")))

#### Plot: VaR by tolerance level (p) ####
data_threshold %>%
  ggplot(aes(x=p,y=Value, color=Method)) +
  geom_line(size=0.8) +
  scale_color_viridis_d(option = "plasma")+
  facet_wrap(~Measure) +
  labs(x="tolerance level (p)",y=NULL,
       title = "VaR by tolerance level (p)") +
  theme_gianluca()



#### Estimate volatility models (EWMA and GARCH), Backtesting and Plot VaR violation ####

ewmafix <- volatmod(y = LogReturn, p = 0.01, v = 3, model = "ewmafix",
                    lambda = 0.94, Tu = 50, date = sp500date)
ewmafix$plot_VaRviol 
ewmafix_test <- ewmafix$backtesting


ewmaig <- volatmod(y = LogReturn, p = 0.01, v = 3, model = "ewmaig",
                   lambda = 0.94, Tu = 50, date = sp500date)
ewmaig$plot_VaRviol
ewmaig_test <- ewmaig$backtesting


garch <- volatmod(y = LogReturn, p = 0.01, v = 3, model = "garch",
                  lambda = 0.94, Tu = 50, date = sp500date)
garch$plot_VaRviol
garch_test <- garch$backtesting


ngarch <- volatmod(y = LogReturn, p = 0.01, v = 3, model = "ngarch",
                   lambda = 0.94, Tu = 50, date = sp500date)
ngarch$plot_VaRviol 
ngarch_test <- ngarch$backtesting



gjrgarch <- volatmod(y = LogReturn, p = 0.01, v = 3, model = "gjr",
                     lambda = 0.94, Tu = 50, date = sp500date)
gjrgarch$plot_VaRviol 
gjrgarch_test <- gjrgarch$backtesting


egarch <- volatmod(y = LogReturn, p = 0.01, v = 3, model = "egarch",
                   lambda = 0.94, Tu = 50, date = sp500date)
egarch$plot_VaRviol 
egarch_test <- egarch$backtesting




#### Select models which pass backtesting ####
bind_rows(ewmafix_test,ewmaig_test,garch_test,
          gjrgarch_test,ngarch_test,egarch_test) %>%
  filter(UCtest=="Fail to Reject H0" & CCtest=="Fail to Reject H0" &
           DURtest=="Fail to Reject H0" & EStest=="Fail to Reject H0") %>%
  select(Model, Method, VaR_Viol=VaR.AViol)


#### Join best models ####
data_best <-  bind_rows(ewmafix$estimates, ewmaig$estimates, garch$estimates,
                        gjrgarch$estimates, ngarch$estimates, egarch$estimates) %>%
  filter(model %in% c("ewmafix","ewmaig","garch","gjr",
                      "ngarch","egarch") & Method %in% c("FHS", "EVT") | 
           model %in% c("gjr","ngarch") & Method == "CF" |
           model == "ewmafix" & Method == "Student") %>%
  mutate_if(is.character,as.factor) %>%
  mutate(model = fct_recode(model, "EWMA (RiskMetrics)"="ewmafix",
                            "EWMA (I-GARCH)"="ewmaig",
                            "GARCH"="garch", "GJR-GARCH"="gjr",
                            "N-GARCH"="ngarch","E-GARCH"="egarch"),
         Method = paste(model,Method,sep = " "))


#### Plot: VaR violation to compare best models ####
data_best %>%
  filter(Measure == "VaR") %>%
  mutate(Violation = ifelse(rt<Value,"VaR violation",
                            "Daily returns SP500")) %>%
  ggplot() +
  geom_point(aes(Day,rt,fill=Violation), shape =21, alpha=0.3) +
  scale_fill_manual(values = c("gray50","red"))+
  geom_line(aes(Day,Value), size = 0.4, color = "cornflowerblue") +
  geom_line(aes(Day,Sigma), color = "chartreuse") +
  facet_wrap(~Method, ncol = 5) + 
  labs(x=NULL, y=NULL,fill=NULL,
       title = "Daily VaR (1%): comparison best models") +
  theme_gianluca()


# Summary 
gjrgarch$fitted_models$empirical

# Plot 
plot(gjrgarch$fitted_models$empirical, which = 8)

# VaR 10 step ahed (Filtered Historical Simulation) 
gjrgarch$bootFHS.10ahed$emp.VaRfhs

# ES 10 step ahed (Filtered Historical Simulation) 
gjrgarch$bootFHS.10ahed$emp.ESfhs


#### Plot: ACF of residuals and std. residuals ####
data.frame(std_res = gjrgarch$fitted_models$empirical@fit$z,
           sq_std_res = gjrgarch$fitted_models$empirical@fit$z^2,
           res = gjrgarch$fitted_models$empirical@fit$residuals) %>%
  pivot_longer(cols = everything(), names_to = "Variable",
               values_to = "Value") %>%
  group_by(Variable) %>%
  summarize(Lag = acf(Value,lag.max = 100,plot = FALSE)$lag[-1],
            ACF = acf(Value,lag.max = 100,plot = FALSE)$acf[-1]) %>%
  ungroup() %>%
  mutate(Variable = fct_relevel(Variable,"sq_std_res",after = Inf),
         Variable = fct_recode(Variable, "Residuals"="res",
                               "Standardized Residuals (z)" = "std_res",
                               "Squared Standardized Residuals (z)" = "sq_std_res")) %>%
  ggplot(aes(x=Lag,y=ACF)) +
  geom_line(color = "cornflowerblue", size = 0.55) +
  facet_wrap(~Variable, ncol = 1, scales = "free") +
  theme_gianluca() +
  labs(title = "Autocorrelations")

  


