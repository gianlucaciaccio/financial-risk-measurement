#### Student Q-Q Plot ####
qqt <- function(y,df=Inf,ylim=range(y),
                main="Student's t Q-Q Plot",
                xlab="Theoretical Quantiles",
                ylab="Sample Quantiles",plot.it=TRUE,...) {
  #	Student's t probability plot
  #	Gordon Smyth
  #	3 Oct 2002
  
  y <- y[!is.na(y)]
  if(0 == (n <- length(y))) stop("y is empty")
  x <- qt(ppoints(n),df=df)[order(order(y))]
  if (plot.it) plot(x,y,main=main,xlab=xlab,ylab=ylab,ylim=ylim,...)
  invisible(list(x=x,y=y))
}





#### THEME GIANLUCA ####
theme_gianluca <-  function(){  
  
  theme(legend.position="bottom", title = element_text(colour = "white"),
        plot.background = element_rect(fill = "gray9",colour = "grey9"),
        legend.background = element_rect(fill = "gray9"),
        legend.key = element_rect(fill = "gray9"),
        legend.text = element_text(colour = "white"),
        legend.title = element_text(colour = "white"),
        strip.background = element_rect(fill = "gray12",colour =  "gray15"),
        panel.background = element_rect(fill = "gray15"), 
        panel.border = element_rect(color = "gray30",fill = NA),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(colour = "gray25"),
        strip.text = element_text(colour = "white"),
        axis.text = element_text(colour = "white"))
}




#### VaR and ES ####
uVaR_ES <- function(y, p = 0.01, type = c("norm", "t", "cf", "evt"),
                    v = NULL,  Tu = NULL, sigma = NULL) {
  
  
  if (length(p) != 1 || is.na(p) || !is.numeric(p) || (p <= 0) || (p > 1)) {
    stop("The argument 'p' must be a single non-NA double value with ", 
         "0 < p < 1.")
  }
  
  
  
  type <- match.arg(type)
  
  if (is.null(sigma)) {
    
    if (length(y) <= 1 || any(is.na(y)) || !is.numeric(y)) {
      stop("A numeric vector of length > 1 and without NAs must be passed to", 
           " 'y'.")
    }
    
    
    s <- sd(y)
    
  } else {
    
    s <- sigma 
    
    if (length(s) <= 1 || any(is.na(s)) || !is.numeric(s) && s <= 0) {
      stop("The argument 'sigma' must be a non-NA double vector with ", 
           "sigma > 0.")
    }
    
  }
  
  
  if (type == "norm") {
    
    VaR   <- s*qnorm(p)          
    ES <- s*-dnorm(qnorm(p))/p
    
  } else if (type == "t") {
    
    if (is.null(v) || length(v) != 1 || is.na(v) || !is.numeric(v)) {
      stop("The argument 'v' must be a single non-NA double value")
    }
    
    Cv <- gamma((v+1)/2)/(gamma(v/2)*sqrt(pi*(v-2)))
    EStp <- (-Cv/p)*((v-2)/(v-1))*(1+(qt(p,v)^2/v))^((1-v)/2)
    
    VaR    <- s*sqrt((v-2)/v)*qt(p,v)
    ES <- s*EStp
    
  } else if (type == "cf") {
    

    sk = skewness(y)
    ek = kurtosis(y)-3
    
    Qn = qnorm(p)
    Qcf = Qn + (sk/6)*(Qn^2-1)+(ek/24)*(Qn^3-3*Qn)-((sk^2)/36)*(2*Qn^3-5*Qn)
    EScf = -((dnorm(Qcf))/p)*(1+(sk/6)*Qcf^3+(ek/24)*(Qcf^4-2*Qcf^2-1))
    
    VaR  <- s*Qcf
    ES <- s*EScf
    
  } else {
    
    if (length(Tu) != 1 || is.na(Tu) || !is.numeric(Tu) || (Tu <= 0)) {
      stop("The argument 'Tu' must be a single non-NA integer value with ", 
           "Tu > 0.")
    }
    
    y <- sort(y, decreasing = TRUE)
    
    Ty <- length(y)
    y_tail <- y[-c(1:(Ty-Tu))]
    u <- min(y[c(1:(Ty-Tu))])
    
    eHill <- sum(log(y_tail/u))/Tu
    
    Fp <- u*(p/(Tu/Ty))^-eHill
    ESevt <- -(u/(eHill-1))*(p/(Tu/Ty))^-eHill
    
    VaR <- as.numeric(s*Fp)
    ES <- as.numeric(s*ESevt)
  }
  
  
  if (is.null(sigma)) {
    
    ES <- mean(y[y<VaR])
    ES <- ifelse(is.nan(ES),VaR,ES)
    
  }
  
  return(list(VaR = VaR, ES = ES))
  
}




#### GARCH VaR ES ####
volatmod <- function(y, date, p, v=NULL, model = c("ewmafix", "ewmaig", "garch",
                                             "ngarch", "gjr", "egarch"),
                     lambda, Tu = NULL) {
  
  model <- match.arg(model)
  

  
  # Set the variance model
  if (model == "ewmafix") {
    
    
    emp.spec <- ugarchspec(variance.model = list(model="iGARCH", garchOrder=c(1,1)),
                           mean.model = list(armaOrder=c(0,0), include.mean=FALSE),
                           fixed.pars = list(alpha1=1-lambda, omega=0))
    
    norm.spec <- ugarchspec(variance.model = list(model="iGARCH", garchOrder=c(1,1)),
                            mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                            fixed.pars = list(alpha1=1-lambda, omega=0), 
                            distribution.model = "norm")
    
    std.spec <- ugarchspec(variance.model = list(model="iGARCH", garchOrder=c(1,1)),
                           mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                           fixed.pars = list(alpha1=1-lambda, omega=0, shape=v), 
                           distribution.model = "std")
    
    main <- "Daily VaR with EWMA (RiskMetrics)"
    
    
  } else if (model == "ewmaig") {
    
    emp.spec <- ugarchspec(variance.model = list(model="iGARCH", garchOrder=c(1,1)),
                           mean.model = list(armaOrder=c(0,0), include.mean=FALSE),
                           fixed.pars = list(omega=0),
                           start.pars = list(alpha1=1-lambda))
    
    norm.spec <- ugarchspec(variance.model = list(model="iGARCH", garchOrder=c(1,1)),
                            mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                            fixed.pars = list(omega=0),
                            start.pars = list(alpha1=1-lambda),
                            distribution.model = "norm")
    
    std.spec <- ugarchspec(variance.model = list(model="iGARCH", garchOrder=c(1,1)),
                           mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                           fixed.pars = list(omega=0, shape=v), 
                           start.pars = list(alpha1=1-lambda),
                           distribution.model = "std")
    
    main <- "Daily VaR with EWMA (I-GARCH)"
    
    
  } else if (model == "garch") {
    
    emp.spec <- ugarchspec(variance.model = list(model="sGARCH", garchOrder=c(1,1)),
                           mean.model = list(armaOrder=c(0,0), include.mean=FALSE))
    
    norm.spec <- ugarchspec(variance.model = list(model="sGARCH", garchOrder=c(1,1)),
                            mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                            distribution.model = "norm")
    
    std.spec <- ugarchspec(variance.model = list(model="sGARCH", garchOrder=c(1,1)),
                           mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                           fixed.pars = list(shape=v), 
                           distribution.model = "std")
    
    main <- "Daily VaR with GARCH"
    
    
  } else if (model == "ngarch") {
    
    
    emp.spec <- ugarchspec(variance.model = list(model="fGARCH", garchOrder=c(1,1), submodel="NAGARCH"),
                           mean.model = list(armaOrder=c(0,0), include.mean=FALSE))
    
    norm.spec <- ugarchspec(variance.model = list(model="fGARCH", garchOrder=c(1,1), submodel="NAGARCH"),
                            mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                            distribution.model = "norm")
    
    std.spec <- ugarchspec(variance.model = list(model="fGARCH", garchOrder=c(1,1), submodel="NAGARCH"),
                           mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                           fixed.pars = list(shape=v), 
                           distribution.model = "std")
    
    
    main <- "Daily VaR with N-GARCH"
    
    
  } else if (model == "gjr") {
    
    emp.spec <- ugarchspec(variance.model = list(model="gjrGARCH", garchOrder=c(1,1)),
                           mean.model = list(armaOrder=c(0,0), include.mean=FALSE))
    
    norm.spec <- ugarchspec(variance.model = list(model="gjrGARCH", garchOrder=c(1,1)),
                            mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                            distribution.model = "norm")
    
    std.spec <- ugarchspec(variance.model = list(model="gjrGARCH", garchOrder=c(1,1)),
                           mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                           fixed.pars = list(shape=v), 
                           distribution.model = "std")
    
    main <- "Daily VaR with GJR-GARCH"
    
    
  } else {
    
    emp.spec <- ugarchspec(variance.model = list(model="eGARCH", garchOrder=c(1,1)),
                           mean.model = list(armaOrder=c(0,0), include.mean=FALSE))
    
    norm.spec <- ugarchspec(variance.model = list(model="eGARCH", garchOrder=c(1,1)),
                            mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                            distribution.model = "norm")
    
    std.spec <- ugarchspec(variance.model = list(model="eGARCH", garchOrder=c(1,1)),
                           mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                           fixed.pars = list(shape=v), 
                           distribution.model = "std")
   
    main <- "Daily VaR with E-GARCH"
    
     
  }
  


  if (model == "ewmafix") {
    
    
    emp.fit <- suppressWarnings(ugarchfit(data = y, spec = emp.spec, solver = 'hybrid'))
    
    norm.fit <- ugarchfit(data = y, spec = norm.spec, solver = 'hybrid')
    
    std.fit <- ugarchfit(data = y, spec = std.spec, solver = 'hybrid')
    
    emp.boot <- ugarchboot(emp.spec, data = y, n.ahead = 10,
                           n.bootpred = 10000, sampling = "raw",
                           rseed = 1:10000, method = "Partial")
    
    norm.boot <- suppressWarnings(ugarchboot(norm.fit, n.ahead = 10, 
                                             n.bootpred = 10000, 
                                             sampling = "raw",
                                             rseed = 1:10000,
                                             method = "Partial"))

    std.boot <- suppressWarnings(ugarchboot(std.fit, n.ahead = 10, 
                                            n.bootpred = 10000,
                                            sampling = "raw",
                                            rseed = 1:10000,
                                            method = "Partial"))
    
    
    emp.z <- as.numeric(emp.fit@filter$z)
    
    
  } else {
    

    emp.fit <- ugarchfit(data = y, spec = emp.spec, solver = 'hybrid')
    
    norm.fit <- ugarchfit(data = y, spec = norm.spec, solver = 'hybrid')
    
    std.fit <- ugarchfit(data = y, spec = std.spec, solver = 'hybrid')
    
    emp.boot <- suppressWarnings(ugarchboot(emp.fit, n.ahead = 10,
                                            n.bootpred = 10000,
                                            sampling = "raw",
                                            rseed = 1:10000,
                                            method = "Partial"))
    
    norm.boot <- suppressWarnings(ugarchboot(norm.fit, n.ahead = 10, 
                                             n.bootpred = 10000, 
                                             sampling = "raw", 
                                             rseed = 1:10000, 
                                             method = "Partial"))
    
    std.boot <- suppressWarnings(ugarchboot(std.fit, n.ahead = 10,
                                            n.bootpred = 10000,
                                            sampling = "raw", 
                                            rseed = 1:10000, 
                                            method = "Partial"))
    
    emp.z <- as.numeric(emp.fit@fit$z)
    
    
  }
  
  
  emp.sigma <- as.numeric(sigma(emp.fit))
  
  norm.sigma <- as.numeric(sigma(norm.fit))
  
  std.sigma <- as.numeric(sigma(std.fit))
  
  
  fitted_models <- list(empirical = emp.fit,
                        normal = norm.fit,
                        student = std.fit)

  
  emp.return10 <- rowSums(emp.boot@fseries)
  
  emp.VaRfhs <- as.numeric(quantile(emp.return10, p))
  
  emp.ESfhs <- sum(emp.return10*(emp.return10 < emp.VaRfhs))/(p*length(emp.return10))

  norm.return10 <- rowSums(norm.boot@fseries)
  
  norm.VaRfhs <- as.numeric(quantile(norm.return10, p))
  
  norm.ESfhs <- sum(emp.return10*(norm.return10 < norm.VaRfhs))/(p*length(norm.return10))
  
  std.return10 <- rowSums(std.boot@fseries)
  
  std.VaRfhs <- as.numeric(quantile(std.return10, p))
  
  std.ESfhs <- sum(emp.return10*(std.return10 < std.VaRfhs))/(p*length(std.return10))
  
  
  bootFHS.10ahed <- list(emp.VaRfhs = emp.VaRfhs, emp.ESfhs = emp.ESfhs,
                         norm.VaRfhs = norm.VaRfhs, norm.ESfhs = norm.ESfhs,
                         std.VaRfhs = std.VaRfhs, std.ESfhs = std.ESfhs)
  
  
  
  Norm_VaR <- uVaR_ES(type = "norm", p = p, sigma = norm.sigma)$VaR
  
  Norm_ES <- uVaR_ES(type = "norm", p = p, sigma = norm.sigma)$ES
  
  Norm_CovTests <- VaRTest(alpha = p, actual = y, VaR = Norm_VaR)
  
  Norm_DurTest <- VaRDurTest(alpha = p, actual = y, VaR = Norm_VaR)
  
  Norm_ESTest <- ESTest(alpha = p, actual = y, ES = Norm_ES, VaR = Norm_VaR)
  
  NormTests <- data.frame(Method="Normal",
                          VaR.EViol=Norm_CovTests$expected.exceed,
                          VaR.AViol=Norm_CovTests$actual.exceed,
                          UCtest=Norm_CovTests$uc.Decision,
                          CCtest=Norm_CovTests$cc.Decision,
                          DURtest=Norm_DurTest$Decision,
                          ES.EViol=Norm_ESTest$expected.exceed,
                          ES.AViol=Norm_ESTest$actual.exceed,
                          EStest=Norm_ESTest$Decision)
  
  Sdt_VaR <- uVaR_ES(type = "t", p = p, sigma = std.sigma, v = v)$VaR
  
  Sdt_ES <- uVaR_ES(type = "t", p = p, sigma = std.sigma, v = v)$ES
  
  Student_CovTests <- VaRTest(alpha = p, actual = y, VaR = Sdt_VaR)
  
  Student_DurTest <- VaRDurTest(alpha = p, actual = y, VaR = Sdt_VaR)
  
  Student_ESTest <- ESTest(alpha = p, actual = y, ES = Sdt_ES, VaR = Sdt_VaR)
  
  StudentTests <- data.frame(Method="Student",
                             VaR.EViol=Student_CovTests$expected.exceed,
                             VaR.AViol=Student_CovTests$actual.exceed,
                             UCtest=Student_CovTests$uc.Decision,
                             CCtest=Student_CovTests$cc.Decision,
                             DURtest=Student_DurTest$Decision,
                             ES.EViol=Student_ESTest$expected.exceed,
                             ES.AViol=Student_ESTest$actual.exceed,
                             EStest=Student_ESTest$Decision)
  
  FHS_VaR <- emp.sigma*quantile(emp.z,p)
  
  FHS_ES <- (emp.sigma*sum(emp.z*(emp.z<FHS_VaR)))/(p*length(emp.z))
  
  FHS_CovTests <- VaRTest(alpha = p, actual = y, VaR = FHS_VaR)
  
  FHS_DurTest <- VaRDurTest(alpha = p, actual = y, VaR = FHS_VaR)
  
  FHS_ESTest <- ESTest(alpha = p, actual = y, ES = FHS_ES, VaR = FHS_VaR)
  
  FHSTests <- data.frame(Method="FHS",
                         VaR.EViol=FHS_CovTests$expected.exceed,
                         VaR.AViol=FHS_CovTests$actual.exceed,
                         UCtest=FHS_CovTests$uc.Decision,
                         CCtest=FHS_CovTests$cc.Decision,
                         DURtest=FHS_DurTest$Decision,
                         ES.EViol=FHS_ESTest$expected.exceed,
                         ES.AViol=FHS_ESTest$actual.exceed,
                         EStest=FHS_ESTest$Decision)
  
  
  CF_VaR <- uVaR_ES(y = emp.z, type = "cf", p = p, sigma = emp.sigma)$VaR
  
  CF_ES <- uVaR_ES(y = emp.z, type = "cf", p = p, sigma = emp.sigma)$ES
  
  CF_CovTests <- VaRTest(alpha = p, actual = y, VaR = CF_VaR)
  
  CF_DurTest <- VaRDurTest(alpha = p, actual = y, VaR = CF_VaR)
  
  CF_ESTest <- ESTest(alpha = p, actual = y, ES = CF_ES, VaR = CF_VaR)
  
  CFTests <- data.frame(Method="CF",
                        VaR.EViol=CF_CovTests$expected.exceed,
                        VaR.AViol=CF_CovTests$actual.exceed,
                        UCtest=CF_CovTests$uc.Decision,
                        CCtest=CF_CovTests$cc.Decision,
                        DURtest=CF_DurTest$Decision,
                        ES.EViol=CF_ESTest$expected.exceed,
                        ES.AViol=CF_ESTest$actual.exceed,
                        EStest=CF_ESTest$Decision)
  
  

  EVT_VaR <- uVaR_ES(y = emp.z, type = "evt", p = p, sigma = emp.sigma, Tu = Tu)$VaR
  
  EVT_ES <- uVaR_ES(y = emp.z, type = "evt", p = p, sigma = emp.sigma, Tu = Tu)$ES
  
  EVT_CovTests <- VaRTest(alpha = p, actual = y, VaR = EVT_VaR)
  
  EVT_DurTest <- VaRDurTest(alpha = p, actual = y, VaR = EVT_VaR)
  
  EVT_ESTest <- ESTest(alpha = p, actual = y, ES = EVT_ES, VaR = EVT_VaR)
  
  EVTTests <- data.frame(Method="EVT",
                         VaR.EViol=EVT_CovTests$expected.exceed,
                         VaR.AViol=EVT_CovTests$actual.exceed,
                         UCtest=EVT_CovTests$uc.Decision,
                         CCtest=EVT_CovTests$cc.Decision,
                         DURtest=EVT_DurTest$Decision,
                         ES.EViol=EVT_ESTest$expected.exceed,
                         ES.AViol=EVT_ESTest$actual.exceed,
                         EStest=EVT_ESTest$Decision)
  
  estimates = list(model = model, 
                   rt = y,
                   zt = emp.z,
                   Empirical_sigma = emp.sigma,
                   FHS_VaR = FHS_VaR,
                   FHS_ES = FHS_ES,
                   CF_VaR = CF_VaR,
                   CF_ES = CF_ES,
                   EVT_VaR = EVT_VaR,
                   EVT_ES = EVT_ES,
                   Normal_sigma = norm.sigma,
                   Normal_VaR = Norm_VaR, 
                   Normal_ES = Norm_ES,
                   Student_sigma = std.sigma,
                   Student_VaR = Sdt_VaR,
                   Student_ES = Sdt_ES)
  
  estimatesDF <- estimates %>%
    bind_cols() %>%
    mutate(Day = date) %>%
    pivot_longer(cols = -c("Day","model","rt","zt","Empirical_sigma",
                           "Normal_sigma","Student_sigma"),
                 names_to = c("Method","Measure"),
                 names_sep = "_",values_to = "Value") %>%
    pivot_longer(cols = c("Empirical_sigma","Normal_sigma",
                          "Student_sigma"),
                 names_to = c("Distribution","NULL"),
                 names_sep = "_", values_to = "Sigma")  %>%
    filter(Method %in% c("CF","FHS","EVT") &  Distribution %in% c("Empirical") |
             Method %in% c("Normal") &  Distribution %in% c("Normal") |
             Method %in% c("Student") & Distribution %in% c("Student")) %>%
    mutate(Method = fct_relevel(Method, c("Normal", "Student", "FHS","CF", "EVT"))) %>%
    select(Day,model,rt,zt,Distribution,Sigma,Method,Measure,Value)
  
  
  
  plot_VaRviol <- estimatesDF %>%
    filter(Measure == "VaR") %>%
    mutate(Violation = ifelse(rt<Value,"VaR violation",
                              "Daily returns SP500")) %>%
    ggplot() +
    geom_point(aes(Day,rt,fill=Violation), shape =21, alpha=0.3) +
    scale_fill_manual(values = c("gray50","red"))+
    geom_line(aes(Day,Value), size = 0.4, color = "cornflowerblue") +
    geom_line(aes(Day,Sigma, color = Distribution)) +
    scale_color_manual(values = c("chartreuse","yellow",
                                  "darkorange")) +
    facet_wrap(~Method, ncol = 3) + 
    labs(y="", title = main) +
    guides(fill = guide_legend(""),
           color = guide_legend(expression("Conditional distribution for "~sigma))) +
    theme_gianluca()
  
  
  
  backtesting <- bind_rows(NormTests, StudentTests, 
                     FHSTests, CFTests, EVTTests) %>%
    mutate(Model = model) %>%
    select(Model, everything())
  
  
  return(list(estimates = estimatesDF, plot_VaRviol = plot_VaRviol,
              backtesting = backtesting, fitted_models = fitted_models,
              bootFHS.10ahed = bootFHS.10ahed))
  
}


