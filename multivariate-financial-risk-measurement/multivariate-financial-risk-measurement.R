library(readxl)
library(zoo)
library(quarks)
library(rugarch)
library(rmgarch)
library(moments)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(cowplot)
library(xts)

source("functions.R")


#### Load the data and calculate returns ####
data <- read_excel("datamulti.xlsx")

data$date <- as.Date(data$date)

sp500ret <- diff(log(data$sp500))

ustnote10ret <- diff(log(data$ustnote10yr))

return_date <- data$date[-1]

assets <- cbind(sp500ret,ustnote10ret)
assets_xts <- as.xts(assets, order.by = return_date)




#### Calculate moments and Jarque-Bera test of log returns ####
data %>%
  mutate(sp500ret = log(sp500/lag(sp500)),
         ustnote10ret = log(ustnote10yr/lag(ustnote10yr))) %>%
  filter(!is.na(sp500ret)) %>%
  select(sp500ret,ustnote10ret) %>%
  pivot_longer(cols = everything(),
               names_to = "Asset",
               values_to = "Value") %>%
  group_by(Asset) %>%
  summarize(mean = mean(Value), sd = sd(Value), 
            skewness = skewness(Value), 
            kurtosis = kurtosis(Value),
            EK = kurtosis(Value)-3)


jarque.test(sp500ret)
jarque.test(ustnote10ret)




#### Plot: histogram log returns and Normal density ####
# SP500
p1 <- data.frame(sp500ret) %>%
  ggplot(mapping = aes(sp500ret)) +
  geom_histogram(aes(y = ..density..), bins = 50,
                 fill = "cornflowerblue", color="gray15") +
  stat_function(fun = dnorm, color="red",
                args = list(mean = mean(sp500ret),
                            sd = sd(sp500ret))) +
  theme_gianluca() +
  labs(x="", y="", title = "SP500") + 
  theme(plot.margin = margin(1,0,0,0, "cm"),
        plot.title = element_text(hjust = 0.5))

# US 10yr T-Note
p2 <- data.frame(ustnote10ret) %>%
  ggplot(mapping = aes(ustnote10ret)) +
  geom_histogram(aes(y = ..density..), bins = 50,
                 fill = "cornflowerblue", color="gray15") +
  stat_function(fun = dnorm, color="red",
                args = list(mean = mean(ustnote10ret),
                            sd = sd(ustnote10ret))) +
  theme_gianluca() +
  labs(x="", y="", title = "US 10yr T-Note") + 
  theme(plot.margin = margin(1,0.2,0,0, "cm"),
        plot.title = element_text(hjust = 0.5))

# Grid
g1 <- plot_grid(p1,p2) + 
  draw_plot_label("Distribution of log returns vs. Normal Density", 
                  color = "white",
                  size = 13,hjust = -0.2)



#### Plot: Q-Q Plot log returns and Normal quantiles ####
# SP500
p3 <- ggplot(mapping = aes(sample=sp500ret)) +
  geom_qq(color="cornflowerblue", shape=1) + 
  geom_qq_line(color="red") +
  theme_gianluca() +
  labs(x="", y="", title = "SP 500") + 
  theme(plot.margin = margin(1,0,0,0, "cm"),
        plot.title = element_text(hjust = 0.5))

# US 10yr T-Note
p4 <- ggplot(mapping = aes(sample=ustnote10ret)) +
  geom_qq(color="cornflowerblue", shape=1) + 
  geom_qq_line(color="red") +
  theme_gianluca() +
  labs(x="", y="", title = "US 10yr T-Note") + 
  theme(plot.margin = margin(1,0.2,0,0, "cm"),
        plot.title = element_text(hjust = 0.5))

# Grid
g2 <- plot_grid(p3,p4) + 
  draw_plot_label("Normal Q-Q Plot of log returns", 
                  color = "white",
                  size = 13,hjust = -0.2) +
  draw_label("Theorical", x=0.5, y=  0, vjust=-0.5,
             angle= 0,color = "white",size = 10) +
  draw_label("Sample", x=  0, y=0.5, vjust= 1.5,
             angle=90,color = "white",size = 10)

#### Grid Histograms and Q-Q Plots
plot_grid(g1,g2,ncol = 1)




#### Plot: Time series of SP500 closing price and log returns ####
# SP500
p5 <- data %>%
  rename(sp500_price = sp500, ustnote10yr_price = ustnote10yr) %>%
  mutate(sp500_return = log(sp500_price/lag(sp500_price)),
         ustnote10yr_return = log(ustnote10yr_price/lag(ustnote10yr_price))) %>%
  pivot_longer(cols = -c("date"), names_to = c("Asset","Index"),
               names_sep = "_", values_to = "Value") %>%
  mutate(Asset = ifelse(Asset=="sp500","SP 500","US 10yr T-Note"),
         Index = ifelse(Index=="price","Closing Price","Log Return")) %>%
  filter(Asset=="SP 500") %>%
  ggplot(aes(x = date, y = Value)) +
  geom_line(color = "cornflowerblue", size = 0.55) +
  facet_grid(Index~Asset,scales = "free") +
  theme_gianluca() +
  theme(strip.background.y = element_blank(),
        strip.text.y =  element_blank()) +
  labs(x = "", y = "", title = "")


# US 10yr T-Note
p6 <- data %>%
  rename(sp500_price = sp500, ustnote10yr_price = ustnote10yr) %>%
  mutate(sp500_return = log(sp500_price/lag(sp500_price)),
         ustnote10yr_return = log(ustnote10yr_price/lag(ustnote10yr_price))) %>%
  pivot_longer(cols = -c("date"), names_to = c("Asset","Index"),
               names_sep = "_", values_to = "Value") %>%
  mutate(Asset = ifelse(Asset=="sp500","SP 500","US 10yr T-Note"),
         Index = ifelse(Index=="price","Closing Price","Log Return")) %>%
  filter(Asset=="US 10yr T-Note") %>%
  ggplot(aes(x = date, y = Value)) +
  geom_line(color = "cornflowerblue", size = 0.55) +
  facet_grid(Index~Asset,scales = "free") +
  theme_gianluca() +
  labs(x = "", y = "", title = "")

# Grid
plot_grid(p5,p6) + 
  draw_plot_label("Closing price and log return",
                  color = "white", size = 12, hjust = -0.1)



#### ACF of log returns and squared log returns ####
data_autocor <- data.frame(
  Lag = acf(sp500ret,lag.max = 100,plot = FALSE)$lag[-1],
  SP500_Return = acf(sp500ret,lag.max = 100,plot = FALSE)$acf[-1],
  SP500_Squared.Return = acf(sp500ret^2,lag.max = 100,plot = FALSE)$acf[-1],
  USTNote10yr_Return = acf(ustnote10ret,lag.max = 100,plot = FALSE)$acf[-1],
  USTNote10yr_Squared.Return = acf(ustnote10ret^2,lag.max = 100,plot = FALSE)$acf[-1]) %>%
  pivot_longer(cols = -c("Lag"), names_to = c("Asset","Return"),
               names_sep = "_", values_to = "ACF") %>%
  mutate(Return = chartr(".", " ", Return),
         Asset = ifelse(Asset=="SP500","SP 500","US 10yr T-Note"))



#### Plot: ACF of log returns and squared log returns ####
# SP500
p7 <- data_autocor %>% filter(Asset=="SP 500") %>%
  ggplot(aes(x=Lag,y=ACF)) +
  geom_line(color = "cornflowerblue", size = 0.55) +
  facet_grid(Return~Asset, scales = "free") +
  theme_gianluca() +
  theme(plot.margin = margin(1,0,0,0, "cm"),
        strip.background.y = element_blank(),
        strip.text.y =  element_blank()) + 
  labs(x = "", y = "")

# US 10yr T-Note
p8 <- data_autocor %>% filter(Asset=="US 10yr T-Note") %>%
  ggplot(aes(x=Lag,y=ACF)) +
  geom_line(color = "cornflowerblue", size = 0.55) +
  facet_grid(Return~Asset, scales = "free") +
  theme_gianluca() + 
  theme(plot.margin = margin(1,0.1,0,0, "cm")) + 
  labs(x = "", y = "")

# Grid
plot_grid(p7,p8) + 
  draw_plot_label("Autocorrelations",color = "white",
                  size = 13,hjust = -0.2) +
  draw_label("Lag", x=0.5, y=  0, vjust=-0.5,
             angle= 0,color = "white",size = 10) +
  draw_label("ACF", x=  0, y=0.5, vjust= 1.5,
             angle=90,color = "white",size = 10)


#### Unconditional Covariance and Correlation ####
cov(sp500ret,ustnote10ret)
cor(sp500ret,ustnote10ret)

#### Static VaR and ES #####
mVaR_ES(y = assets, p = 0.01, type = "hs", w = c(0.5,0.5))
mVaR_ES(y = assets, p = 0.01, type = "whs", lambda = 0.94, w = c(0.5,0.5))
mVaR_ES(y = assets, p = 0.01, type = "norm", w = c(0.5,0.5))
mVaR_ES(y = assets, p = 0.01, type = "t", v=3, w = c(0.5,0.5))
mVaR_ES(y = assets, p = 0.01, type = "cf", w = c(0.5,0.5))
mVaR_ES(y = assets, p = 0.01, type = "evt", Tu = 50, w = c(0.5,0.5))


#### Test of Dynamic Correlation ####
DCCtest(assets)


#### DCC-Garch (Standard, NGARCH and EWMA-RiskMetrics), VaR and ES ####

#### EWMA ####
ewma.dcc <- TdccVaR_ES(y = assets, date = return_date, w = c(0.5,0.5),
                       p = 0.01, v = 3, model = "ewma", lambda = 0.94)

ewma.dcc$plotDCC
ewma.dcc$Tests

#### GARCH ####
garch.dcc <- TdccVaR_ES(y = assets, date = return_date, w = c(0.5,0.5),
                        p = 0.01, v = 3, model = "garch")

garch.dcc$plotDCC
garch.dcc$Tests
garch.dcc$model.fit

#### NGARCH ####
ngarch.dcc <- TdccVaR_ES(y = assets, date = return_date, w = c(0.5,0.5),
                         p = 0.01, v = 3, model = "ngarch")

ngarch.dcc$plotDCC
ngarch.dcc$Tests
ngarch.dcc$model.fit


