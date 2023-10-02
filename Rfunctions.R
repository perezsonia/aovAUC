library(MASS)

# Function: vAUC ----

# Objective: AUC estimation (hat{A}) and its variance (hat{sigma}^2)
# Arguments: 
## V = (bio)marker values
## D = population indicator values where 0 = negative (control) and 1 = positive (case)
# Output: List of 2
## AUC = AUC estimation
## vAUC = sigma^2 estimation (in page 5)

vAUC <- function(V, D){
  
 V0 <- V[D==0]; V1 <- V[D==1]
 n0 <- length(V0); n1 <- length(V1)
 F0 <- ecdf(V0)(V1); F1 <- ecdf(V1)(V0)
 S01 <- max(var(F0), 1e-3); S10 <- max(var(F1), 1e-3)
 
 return(list(AUC = mean(F0), vAUC = S10/n0 + S01/n1))
 
}


# Function: reV ----

# Objective: Inter-group variance estimation, hat{tau}^2 (Appendix)
# Arguments: 
## auc = AUC values (output of vAUC function)
## var = AUC variance (output of vAUC function)
# Output: List of 2
## tau = inter-group variance estimation
## SD = standard deviation estimation of AUC values

reV <- function(auc, var){
  
 V <- var(auc)
 tau <- max(V - mean(var), 0)
 list(tau = tau, SD = V^0.5) 
 
}

# Function: subjects_data ----

# Objective: AUC estimation and variance related measures per subject
# Argument: data = data frame with 4 variables: 
## ID = subject ID
## TrT = treatment
## Gr = population indicator values where 0 = negative (control) and 1 = positive (case)
## Values = (bio)marker values
# Output: Data frame of 8 variables where each row represents each subject
## TrT = treatment
## n0 = number of pre-treatment (control) measures, n1 = number of post-treatment (case) measures
## AUC = AUC estimation (vAUC function output)
## vAUC = sigma^2 estimation (vAUC function output)
## vE = inter-group variance estimation (reV function output)
## SD = standard deviation estimation of AUC values (reV function output)
## vToT = total variance estimation of AUC = vAUC + vE

subjects_data <- function(data){

   rats <- as.data.frame(
      do.call(rbind, by(data, data$ID, function(x){
         TrT <- mean(x$TrT)
         A <- vAUC(x$Values, x$Gr)
         n0 <- sum(x$Gr == 0); n1 <- sum(x$Gr == 1)
         c(TrT = TrT, n0 = n0, n1 = n1, AUC = A$AUC, vAUC = A$vAUC)
         }))
      )

   tau.TrT <- do.call(rbind, by(rats, rats$TrT, function(x){
      inter <- reV(x$AUC, x$vAUC)
      c(vE = inter$tau, SD = inter$SD)
   }))
   rats <- cbind(rats, tau.TrT[rats$TrT,])
   rats$vToT <- rats$vAUC + rats$vE
   
   return(rats)
   
}



# Function: summary.subjectsdata ----

# Objective: Summary for output of subjects_data function
# Argument: subjects.data = data frame output of subjects_data function

summary.subjectsdata <- function(subjects.data){
  
  Table <- do.call(rbind, by(subjects.data, subjects.data$TrT, function(x){
    mean <- mean(x$AUC)
    sd <- unique(x$SD)
    min <- min(x$AUC); max <- max(x$AUC)
    n <- nrow(x)
    mP <- mean(x$n1); mN <- mean(x$n0)
    c(round(c(mean, sd, min, max), 2), n, round(c(mP, mN)))
  }))
  
  rownames(Table) <- paste("Treatment", levels(as.factor(subjects.data$TrT)))
  colnames(Table) <- c(paste0(c("Mean", "Sd", "Min", "Max"), ".AUC"), c("n", "m.mP", "m.mN"))
  Table
  
}


# Function: SSE ----

# Objective: Sum of the intra-groups variability
# Argument: rats = data frame output of subjects_data function
# Output: List of 3
## summands = summands for each treatment/group used to calculate the sum of the intra-groups variability (eq. (5))
## SSE = sum of the intra-groups variability
## df = degrees of freedom for SSE distribution (chi-squared)

SSE <- function(rats){
   
   treatments <- unique(rats$TrT)
   
   SSEk <- sapply(treatments,  function(trt){
      x <- subset(rats, rats$TrT == trt)
      nd <- nrow(x)
      V <- diag(x$vToT)
      M <- diag(rep(1,nd)) - matrix(1/nd, ncol = nd, nrow = nd)
      P <- ginv(M %*% V %*% M) # Moore-Penrose generalized inverse of matrix M%*%V%*%M
      t(x$AUC) %*% P %*% x$AUC
   })
   
   list(summands = SSEk, SSE = sum(SSEk), df = nrow(rats) - length(treatments))
}


# Function: SSF ----

# Objective: Sum of the inter-groups variability
# Argument: rats = data frame output of subjects_data function
# Output: List of 2
## SSF = sum of the inter-groups variability
## df = degrees of freedom for SSF distribution (chi-squared)

SSF <- function(rats){
   
   k <- length(table(rats$TrT))
   trtMeans <- tapply(rats$AUC, rats$TrT, mean)
   eps <- (mean(rats$AUC)/0.6)^0.25
   trtVars <- tapply(rats$vToT, rats$TrT, mean) / table(rats$TrT) * eps
   V <- diag(trtVars)
   M <- diag(rep(1,k)) - matrix(1/k, ncol = k, nrow = k)
   P <- ginv(M %*% V %*% M) # Moore-Penrose generalized inverse of matrix M%*%V%*%M
   
   list(SSF = as.numeric(t(trtMeans) %*% P %*% trtMeans), df = length(trtVars) - 1)
   
}


# Function: aucAOV ----

# Objective: ANOVA-type hypothesis testing for comparing the difference in the impact between several treatments
# Argument: rats = data frame output of subjects_data function
# Output: List of 5
## Means = AUC estimation mean for each treatment
## F.value = test statistic value
## p.value = resulting p-value
## SumOfSquares = list of 4: SSE, SSE.df, SSF and SSF.df (output of SSE and SSF functions)
## iAUC = input data frame

aucAOV <- function(rats){
   
   trtMeans <- tapply(rats$AUC, rats$TrT, mean)
   sse <- SSE(rats); ssf <- SSF(rats)
   F.value <- (ssf$SSF / ssf$df) / (sse$SSE / sse$df)
   p.value <-  1 - pf(F.value, ssf$df, sse$df)
   
   list(Means = trtMeans, F.value = F.value, p.value = p.value,
        SumOfSquares = c(SSE = sse$SSE, SSE.df = sse$df, SSF = ssf$SSF, SSF.df = ssf$df),
        iAUC = rats)
   
}


# Function: HST ----

# Objective: Distribution approximation for Algorithm Tukey's HSD-type test (post hoc) statistic
# Arguments: 
## k = number of treatments to compare
## B = number of generated values for post hoc statistic mD[k]
# Output: Vector of mD generated values (length B)

HST <- function(k, B = 1e7){
  
   npairs <- k*(k-1)/2
   X <- matrix(rnorm(k*B), ncol = k, nrow = B)
   Dist <- matrix(NA, ncol = npairs, nrow = B) 
   l <- 1
   for (j in 1:(k-1)) {
      for (i in (j+1):k) {
         Dist[,l] <- abs(X[,j] - X[,i])/sqrt(2)
         l <- l+1}}
   
   apply(Dist, 1, max)
   
}



# Function: posthoc ----

# Objective: Tukey HSD-type hypothesis testing for comparing the difference in the impact between all pairs of treatments
# Arguments: 
## rats = data frame output of subjects_data function
## alpha = significance level
## DIST = vector output of HST(k = length(treatments), B = B)
## ecdf.DIST = empirical cumulative distribution function of DIST = ecdf(DIST)
## B = number of generated values for post hoc statistic mD[k] (only used if DIST and ecdf.DIST are not provided)
# Output: List of 5
## difference = difference between the pairs (Delta[k] in page 8)
## p.values = p values
## threshold = threshold used for identifying statistical significance
## rejection = proportion of rejections
## treatments = treatments

posthoc <- function(rats, alpha = .05, DIST = NULL, B = 1e5, ecdf.DIST = NULL){
   
   nrat <- as.numeric(table(rats$TrT))

   treatments <- as.numeric(levels(as.factor(rats$TrT)))
   pairs <- combn(treatments,2)
   npairs <- ncol(pairs)
   names <- apply(pairs, 2, function(x) paste0(x[1],"-",x[2]))
   
   trtMeans <- tapply(rats$AUC, rats$TrT, mean)
   trtVars <- tapply(rats$vToT, rats$TrT, mean)
   
   difference <- sapply(1:npairs, function(pair.index){
      i <- pairs[1,pair.index]; j <- pairs[2,pair.index]
      abs(trtMeans[j] - trtMeans[i]) / sqrt( trtVars[j]/nrat[j] + trtVars[i]/nrat[i])
      
   })
   names(difference) <- names
   
   if(is.null(ecdf.DIST)){
      DIST <- HST(k = length(treatments), B = B)
      ecdf.DIST <- ecdf(DIST)
   }
   
   p.values <- sapply(1:npairs, function(pair.index) 1 - ecdf.DIST(difference[pair.index]))
   names(p.values) <- names
   
   threshold <- quantile(DIST, 1 - alpha)
   
   list(difference = difference, p.values = p.values, threshold = threshold, 
        rejection = (difference > threshold), treatments = treatments)

}



# Function: aovAUC ----

# Objective: ANOVA-type hypothesis testing for comparing the difference in the impact between several treatments.
#             It also includes Tukey HSD-type hypothesis testing for providing full pair comparisons.

# Arguments:
## formula: Marker ~ Treatment
## id.var = `variable name modelling subject ID`
## group.var = `variable name modelling negative vs. positive`
## data = full data frame
## ph = logical value indicating if the post hoc Tukey's HSD-type test should be performed
## alpha = significance level
## DIST = vector output of HST(k = length(treatments), B = B)
## ecdf.DIST = empirical cumulative distribution function of DIST = ecdf(DIST)
## B = number of generated values for post hoc statistic mD[k] (only used if DIST and ecdf.DIST are not provided)

# Output: List of 12 (object of class `aovauc`)
## Means = AUC estimation mean for each treatment
## F.value = test statistic value
## p.value = resulting p-value
## SumOfSquares = list of 4: SSE, SSE.df, SSF and SSF.df (output of SSE and SSF functions)
## iAUC = input data frame
## RandomEffectsSE = inter-group standard deviation estimation per treatment
## difference = difference between the pairs (Delta[k] in page 8)
## p.values = p values
## threshold = threshold used for identifying statistical significance
## rejection = proportion of rejections
## treatments = treatments
## call = input formula

aovAUC <- function(formula, id.var = "ID", group.var = "Gr", data, ph = TRUE, alpha = .05, DIST = NULL, ecdf.DIST = NULL, B = 1e5){
   
   dt <- data.frame(ID = subset(data, select = as.character(id.var)),
                    Gr = subset(data, select = as.character(group.var)),
                    Values = subset(data, select = as.character(as.list(formula)[[2]])),
                    TrT = subset(data, select = as.character(as.list(formula)[[3]])))
   
   dt.auc <- subjects_data(dt)
   
   result <- aucAOV(dt.auc)
   
   result$RandomEffectsSE <- round(sqrt(tapply(dt.auc$vE, dt.auc$TrT, unique)), 3)
   
   if(ph){
      
      result <- c(result, posthoc(dt.auc, alpha = alpha, DIST = DIST, ecdf.DIST = ecdf.DIST, B = B))
      
   }
   
   result$call <- formula
   
   class(result) <- "aovauc"
   
   result
   
}



# Function: print.aovauc ----

# Objective: Print for aovAUC function
# Argument: x = object of class `aovauc` output of aovAUC function

print.aovauc <- function(x){
   
   cat("\nCall:\n", paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep = "")
   
   Summ1.temp <- cbind(`Sum Square` = round(x$SumOfSquares[c("SSE", "SSF")],2),
                       DF = x$SumOfSquares[c("SSE.df", "SSF.df")])
   Summ1 <- cbind(Summ1.temp,
                  `Mean Square` = round(Summ1.temp[,1]/Summ1.temp[,2],2),
                  `F-Snedecor` = c(round(x$F.value,4),NA),
                  `p-value` = c(x$p.value,NA))
   Summ1 <- rbind(Summ1, c(sum(Summ1.temp[,1]), rep(NA,4)))
   rownames(Summ1) <- c("Intra-group", "Inter-groups", "Total")
   printCoefmat(Summ1, na.print = "", cs.ind = 2)
   
   cat("\nAverage random-effects standard error of ", round(mean(x$RandomEffectsSE),3), " (", paste0(x$RandomEffectsSE, collapse = ", "),")\n\n", sep = "")
   
   cat(paste(rep("-",70), collapse= ""))
   
   cat("\n\nPost hoc test (p-values)\n\n")
   
   PH.comparisons <- sort(names(x$p.values))
   PH.p.values <- x$p.values[order(names(x$p.values))]
   posit.PH.comparison <- lapply(strsplit(PH.comparisons, split = "-"), as.numeric)
   trts <- x$treatments
   Summ2 <- matrix(NA, nrow = length(trts), ncol = length(trts), dimnames = list(paste("Treatment", trts), paste("Treatment", trts)))
   for(comp in seq_along(posit.PH.comparison)){
      i <- posit.PH.comparison[[comp]][1]; j <- posit.PH.comparison[[comp]][2]
      if(i < j) Summ2[i, j] <- PH.p.values[comp] else Summ2[j, i] <- PH.p.values[comp]
   }
   Summ2 <- Summ2[-nrow(Summ2),]
   printCoefmat(round(Summ2,4), na.print = "", has.Pvalue = TRUE)
   
}

