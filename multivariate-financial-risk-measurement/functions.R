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
mVaR_ES <- function(y, p = 0.01, type = c("hs", "whs", "norm", "t", "cf", "evt"),
                    lambda = NULL, v = NULL,  Tu = NULL, w = c(0.5,0.5),
                    sigmaPF = NULL) {
  
  
  
  
  
  type <- match.arg(type)
  
  if (length(p) != 1 || is.na(p) || !is.numeric(p) || (p <= 0) || (p > 1)) {
    stop("The argument 'p' must be a single non-NA double value with ", 
         "0 < p < 1.")
  }
  
  
  if (length(w) <= 1 || any(is.na(w)) || !is.numeric(w)) {
    stop("A numeric vector of length > 1 and without NAs must be passed to", 
         " 'w'.")
  }
  
  if (sum(w) != 1 ) {
    stop("The sum of weights must be equal to 1")
  }
  
  
  
  if (is.null(sigmaPF)) {
    
    if (!is.matrix(y) || length(y) <= 1 || any(is.na(y)) || !is.numeric(y)) {
      stop("A numeric matrix of length > 1 and without NAs must be passed to", 
           " 'y'.")
    }
    
    if (length(w) != ncol(y)) {
      stop("Number of weights must be equal to the number of assets (ncol(y))")
    }
    
    
    sPF <- sqrt(t(w)%*%cov(y)%*%w)
    
    
    if (type == "hs" || type == "whs" || type == "cf" || type == "evt") {
      
      y_PF <- rowSums(w*y)
      
    }
    
    
  } else {
    
    
    if (type == "hs" || type == "whs" || type == "cf" || type == "evt") {
      
      if (!is.matrix(y) || length(y) <= 1 || any(is.na(y)) || !is.numeric(y)) {
        stop("A numeric matrix of length > 1 and without NAs must be passed to", 
             " 'y'.")
      }
      
      y_PF <- rowSums(w*y)
      
    }
    
    sPF <- as.numeric(sigmaPF) 
    
    if (length(sPF) <= 1 || any(is.na(sPF)) || !is.numeric(sPF)) {
      stop("The argument 'sigmaPF' must be a non-NA double vector with", 
           "sigma > 0.")
    }
    
  }
  
  
  if (type == "hs") {
    
    VaR <- -hs(y_PF,p = 1-p, method = "plain")$VaR
    ES <- -hs(y_PF,p = 1-p, method = "plain")$ES
    
  } else if (type == "whs") {
    
    if (length(lambda) != 1 || is.na(lambda) || !is.numeric(lambda) || (lambda <= 0) || (lambda > 1)) {
      stop("The argument 'lambda' must be a single non-NA double value with ", 
           "0 < lambda < 1.")
    }
    
    
    VaR <- -hs(y_PF,p = 1-p, method = "age", lambda = lambda)$VaR
    ES <- -hs(y_PF,p = 1-p, method = "age", lambda = lambda)$ES
    
  } else if (type == "norm") {
    
    VaR   <- as.numeric(sPF*qnorm(p))          
    ES    <- as.numeric(sPF*-dnorm(qnorm(p))/p) 
    
    
  } else if (type == "t") {
    
    if (is.null(v) || length(v) != 1 || is.na(v) || !is.numeric(v)) {
      stop("The argument 'v' must be a single non-NA double value")
    }
    
    Cv <- gamma((v+1)/2)/(gamma(v/2)*sqrt(pi*(v-2)))
    EStp <- (-Cv/p)*((v-2)/(v-1))*(1+(qt(p,v)^2/v))^((1-v)/2)
    
    VaR    <- as.numeric(sPF*sqrt((v-2)/v)*qt(p,v))
    ES     <- as.numeric(sPF*EStp)
    
  } else if (type == "cf") {
    
    sk = skewness(y_PF)
    ek = kurtosis(y_PF)-3
    
    Qn = qnorm(p)
    Qcf = Qn + (sk/6)*(Qn^2-1)+(ek/24)*(Qn^3-3*Qn)-((sk^2)/36)*(2*Qn^3-5*Qn)
    EScf = -((dnorm(Qcf))/p)*(1+(sk/6)*Qcf^3+(ek/24)*(Qcf^4-2*Qcf^2-1))
    
    VaR  <- as.numeric(sPF*Qcf)
    ES <- as.numeric(sPF*EScf)
    
    if (is.null(sigmaPF)) {
     
      ES <- mean(y[y<VaR]) 
    }
    
  } else {
    
    if (length(Tu) != 1 || is.na(Tu) || !is.numeric(Tu) || (Tu <= 0)) {
      stop("The argument 'Tu' must be a single non-NA integer value with ", 
           "Tu > 0.")
    }
    
    
    y_PF <- sort(y_PF, decreasing = TRUE)
    
    Ty <- length(y_PF)
    y_tail <- y_PF[-c(1:(Ty-Tu))]
    u <- min(y_PF[c(1:(Ty-Tu))])
    
    eHill <- sum(log(y_tail/u))/Tu
    
    Fp <- u*(p/(Tu/Ty))^-eHill
    ESevt <- -(u/(eHill-1))*(p/(Tu/Ty))^-eHill
    
    VaR <- as.numeric(sPF*Fp)
    ES <- as.numeric(sPF*ESevt)
    
    
  }
  
  
  return(list(VaR = VaR, ES = ES))
  
}






#### DCC VaR and ES with t distribution ####
TdccVaR_ES <- function(y, date, w = c(0.5,0.5), p, v,
                       model = c("ewma", "garch", "ngarch"), lambda) {
  
  
  
  if (length(p) != 1 || is.na(p) || !is.numeric(p) || (p <= 0) || (p > 1)) {
    stop("The argument 'p' must be a single non-NA double value with ", 
         "0 < p < 1.")
  }
  
  
  model <- match.arg(model)
  
  
  if (length(w) <= 1 || any(is.na(w)) || !is.numeric(w)) {
    stop("A numeric vector of length > 1 and without NAs must be passed to", 
         " 'w'.")
  }
  
  if (sum(w) != 1 ) {
    stop("The sum of weights must be equal to 1")
  }
  
  if (length(w) != ncol(y)) {
    stop("Number of weights must be equal to the number of assets (ncol(y))")
  }
  
  
  if (is.null(v) || length(v) != 1 || is.na(v) || !is.numeric(v)) {
    stop("The argument 'v' must be a single non-NA double value")
  }
  
  
  
  
  if (model == "ewma") {
    
    if (length(lambda) != 1 || is.na(lambda) || !is.numeric(lambda) || 
        (lambda <= 0) || (lambda > 1)) {
      stolambda("The argument 'lambda' must be a single non-NA double value with ", 
                "0 < lambda < 1.")
    }
    
    
    
    main = "DCC-EWMA"
    
    uspec <- ugarchspec(variance.model = list(model="iGARCH", garchOrder=c(1,1)),
                        mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                        fixed.pars = list(alpha1=1-lambda, omega=0, shape=v), 
                        distribution.model = "std")
    
  } else if (model == "garch") {
    
    main = "DCC-GARCH"
    
    uspec <- ugarchspec(variance.model = list(model="sGARCH", garchOrder=c(1,1)),
                        mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                        fixed.pars = list(shape=v), distribution.model = "std")
    
    
  } else {
    
    main = "DCC-NGARCH"
    
    uspec <- ugarchspec(variance.model = list(model="fGARCH", garchOrder=c(1,1), submodel="NAGARCH"),
                        mean.model = list(armaOrder=c(0,0), include.mean=TRUE),
                        fixed.pars = list(shape=v), distribution.model = "std")
    
  }
  
  y_xts <- as.xts(assets, order.by = date)
  
  mspec <- multispec(replicate(2, uspec))
  dccspec <- dccspec(mspec, dccOrder = c(1, 1),
                     model = "DCC",distribution = 'mvt')
  
  mfit <- multifit(mspec, data = y_xts)
  suppressWarnings(dccfit <- dccfit(dccspec, data = y_xts, fit = mfit))
  
  sigma_DCC <- sigma(dccfit)
  colnames(sigma_DCC) <- paste(colnames(sigma_DCC),"sigma", sep = "_")
  
  cov_DCC <- rcov(dccfit, output="matrix")[,2]
  colnames(cov_DCC) <- "covariance"
  
  cor_DCC <- rcor(dccfit, output="matrix")
  colnames(cor_DCC) <- "correlation"
  
  sigmaPF_DCC <- rowSums(w^2*sigma_DCC)+2*prod(w)*cor_DCC*apply(sigma_DCC, 1, prod)
  colnames(sigmaPF_DCC) <- "sigmaPF"
  
  data_DCC <- fortify.zoo(merge.xts(y_xts, cov_DCC, cor_DCC),
                          names = c(Index = "date")) %>%
    mutate(returnPF = rowSums(w*y),
           VaR = mVaR_ES(w = w, p = p, type = "t", v = v,
                         sigmaPF = sigmaPF_DCC)$VaR,
           ES = mVaR_ES(w = w, p = p, type = "t", v = v,
                        sigmaPF = sigmaPF_DCC)$ES) %>%
    select(date,returnPF,covariance,correlation,VaR,ES)
  
  
  p1 <- data_DCC %>%
    select(date,covariance,correlation) %>%
    pivot_longer(cols = -c(date),names_to = "measure") %>%
    mutate(measure = paste("Conditional",measure,sep = " ")) %>%
    ggplot(aes(x=date,y=value)) +
    geom_line(color = "cornflowerblue", size = 0.55) +
    facet_wrap(~measure,scales = "free", ncol = 1) +
    theme_gianluca() +
    theme(plot.margin = margin(1,0,0,0, "cm")) + 
    labs(x = "", y = "")
  
  
  p2 <- data_DCC %>%
    select(date,returnPF,VaR,ES) %>%
    mutate(violation = ifelse(returnPF<VaR,"VaR violations","PF return")) %>%
    pivot_longer(cols = -c(date,returnPF,violation),
                 names_to = "measure") %>%
    ggplot() +
    geom_point(aes(x = date, y = returnPF, fill = violation), 
               shape = 21, alpha = 0.2) +
    scale_fill_manual(values = c("gray50","red")) +
    geom_line(aes(x = date, y = value, color = measure), size = 0.55) +
    scale_color_manual(values = c("cornflowerblue","chartreuse")) +
    theme_gianluca() +
    theme(plot.margin = margin(1,0.2,0,0, "cm")) + 
    labs(x = "", y = "", title = "VaR violations",
         color = "", fill = "")
  
  
  plotDCC <- plot_grid(p1, p2, ncol = 2) + 
    draw_plot_label(main, color = "white", size = 15, hjust = -0.5)
  
  
  test1 <- VaRTest(alpha = 0.01,actual = data_DCC$returnPF, VaR = data_DCC$VaR)
  test2 <- VaRDurTest(alpha = 0.01,actual = data_DCC$returnPF, VaR = data_DCC$VaR)
  test3 <- ESTest(alpha = 0.01,actual = data_DCC$returnPF, ES = data_DCC$ES,
                  VaR = data_DCC$VaR, boot = FALSE)
  
  
  
  Tests <- data.frame(Model=model,
                      VaR.ee=test1$expected.exceed,
                      VaR.ae=test1$actual.exceed,
                      UCtest=test1$uc.Decision,
                      CCtest=test1$cc.Decision,
                      DURtest=test2$Decision,
                      ES.ee=test3$expected.exceed,
                      ES.ae=test3$actual.exceed,
                      EStest=test3$Decision)
  
  
  
  return(list(model.fit =dccfit, plotDCC = plotDCC, Tests = Tests))
  
}




