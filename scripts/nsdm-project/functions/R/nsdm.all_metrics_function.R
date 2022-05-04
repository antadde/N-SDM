#' Evaluation function (define some confusion matrix indices)
#'
#' @author Philipp Brun (philipp.brun@wsl.ch) and Antoine Adde (aadde@unil.ch)

# define some confusion matrix indices

# Sensitivity
sens=function(a,c){a/(a+c)}

# Specificity
spec=function(b,d){d/(b+d)}

# Positive predictive value
ppv=function(a,b){a/(a+b)}

# Negative predictive value
npv=function(c,d){d/(c+d)}

# Jaccard index
jaccard=function(a,b,c){a/(a+b+c)}

# True skill statistic
tss=function(a,b,c,d){a/(a+c)+d/(b+d)-1}

# Accuracy
acc=function(a,b,c,d){(a+d)/(a+b+c+d)}

# Approximations
minusInf <- function(a = a, b = b, c = c, d = d) {
  list2env(list(a = b, b = a, c = d, d = c), envir = .GlobalEnv)
  return(-1e+09)
}

plusZero <- function(a = a, b = b, c = c, d = d) {
  list2env(list(a = b, b = a, c = d, d = c), envir = .GlobalEnv)
  return(1e-09)
}

# helper
# prevalence or base rate
BaseRate <- function(a = a, b = b, c = c, d = d) {
list2env(list(a = b, b = a, c = d, d = c), envir = .GlobalEnv)
  return ((a + c) / (a + b + c + d))
}

# helper
# sensitivity or hit rate
HitRate <- function(a = a, b = b, c = c, d = d) {
list2env(list(a = b, b = a, c = d, d = c), envir = .GlobalEnv)
  return (a / (a + c))
}

# 1 minus hit rate
OneMHR <- function(a = a, b = b, c = c, d = d) {
list2env(list(a = b, b = a, c = d, d = c), envir = .GlobalEnv)
  if (isTRUE(HitRate(a = a, b = b, c = c, d = d) == 1)) {
    return(plusZero(a = a, b = b, c = c, d = d))
  }
  else
    return(1 - (a / (a + c)))
}

# log of hit rate
logHR <- function(a = a, b = b, c = c, d = d) {
list2env(list(a = b, b = a, c = d, d = c), envir = .GlobalEnv)
  if (isTRUE(HitRate(a = a, b = b, c = c, d = d) == 0)) {
    return(minusInf(a = a, b = b, c = c, d = d))
  }
  else
    return(log(HitRate(a = a, b = b, c = c, d = d)))
}

# helper
# false alarm
FalseAlarm <- function(a = a, b = b, c = c, d = d) {
list2env(list(a = b, b = a, c = d, d = c), envir = .GlobalEnv)
  return (b / (b + d))
}

# 1 minus false alarm
OneMFA <- function(a = a, b = b, c = c, d = d) {
list2env(list(a = b, b = a, c = d, d = c), envir = .GlobalEnv)
  if (isTRUE(FalseAlarm(a = a, b = b, c = c, d = d) == 1)) {
    return(plusZero(a = a, b = b, c = c, d = d))
  }
  else
    return (1 - (b / (b + d)))
}

# log of false alarm
logFA <- function(a = a, b = b, c = c, d = d) {
list2env(list(a = b, b = a, c = d, d = c), envir = .GlobalEnv)
  if (isTRUE(FalseAlarm(a = a, b = b, c = c, d = d) == 0)) {
    return(minusInf(a = a, b = b, c = c, d = d))
  }
  else
    return(log(FalseAlarm(a = a, b = b, c = c, d = d)))
}

# main function
# symmetric extremal dependence index
sedi <- function(a = a, b = b, c = c, d = d) {
list2env(list(a = b, b = a, c = d, d = c), envir = .GlobalEnv)

  if (isTRUE(x = (max(a,b,c,d) <= 1) && sum(a,b,c,d) <= 1)) { # no value can be computed
    return(list(NaN))  
  }
  
  if (isTRUE(x = (min(a,b,c,d) < 1))) { # catch all matrices which contain at least one zero
    if (isTRUE(x = (sum(a,b) < 1) &&  (sum(c,d) > 1) )) {
      return(list(0))
    }
    if (isTRUE(x = (sum(c,d) < 1) &&  (sum(a,b) > 1) )) {
      return(list(0))
    }
    if (isTRUE(x = (sum(b,c) < 1) &&  (sum(a,d) > 1) )) {
      return(list(1))
    }
    if (isTRUE(x = (sum(a,d) < 1) &&  (sum(b,c) > 1) )) {
      return(list(0))
    }
  
    if (isTRUE(x = (a < 1 && min(b,c,d) >=1))) { # only true presences equal to zero
      return(list(
                  (logFA(a = a+1e-09, b = b, c = c, d = d) - logHR(a = a+1e-09, b = b, c = c, d = d) - log(OneMFA(a = a+1e-09, b = b, c = c, d = d)) + log(OneMHR(a = a+1e-09, b = b, c = c, d = d)))
                  /
                  (logFA(a = a+1e-09, b = b, c = c, d = d) + logHR(a = a+1e-09, b = b, c = c, d = d) + log(OneMFA(a = a+1e-09, b = b, c = c, d = d)) + log(OneMHR(a = a+1e-09, b = b, c = c, d = d)))))  
    }
  
    if (isTRUE(x = (d < 1 && min(a,b,c) >=1))) { # only true absences (also pseudo-absences, incl. background) equal to zero
      return(list(
                  (logFA(a = a, b = b, c = c, d = d+1e-09) - logHR(a = a, b = b, c = c, d = d+1e-09) - log(OneMFA(a = a, b = b, c = c, d = d+1e-09)) + log(OneMHR(a = a, b = b, c = c, d = d+1e-09)))
                  /
                  (logFA(a = a, b = b, c = c, d = d+1e-09) + logHR(a = a, b = b, c = c, d = d+1e-09) + log(OneMFA(a = a, b = b, c = c, d = d+1e-09)) + log(OneMHR(a = a, b = b, c = c, d = d+1e-09)))))  
    }
  
    if (isTRUE(x = (b < 1 && min(a,c,d) >=1))) { # only zero commission errors
      return(list(
                  (logFA(a = a, b = b+1e-09, c = c, d = d) - logHR(a = a, b = b+1e-09, c = c, d = d) - log(OneMFA(a = a, b = b+1e-09, c = c, d = d)) + log(OneMHR(a = a, b = b+1e-09, c = c, d = d)))
                  /
                  (logFA(a = a, b = b+1e-09, c = c, d = d) + logHR(a = a, b = b+1e-09, c = c, d = d) + log(OneMFA(a = a, b = b+1e-09, c = c, d = d)) + log(OneMHR(a = a, b = b+1e-09, c = c, d = d)))))  
    }
  
    if (isTRUE(x = (c < 1 && min(a,b,d) >=1))) { # only zero omission errors
      return(list(
                  (logFA(a = a, b = b, c = c+1e-09, d = d) - logHR(a = a, b = b, c = c+1e-09, d = d) - log(OneMFA(a = a, b = b, c = c+1e-09, d = d)) + log(OneMHR(a = a, b = b, c = c+1e-09, d = d)))
                  /
                  (logFA(a = a, b = b, c = c+1e-09, d = d) + logHR(a = a, b = b, c = c+1e-09, d = d) + log(OneMFA(a = a, b = b, c = c+1e-09, d = d)) + log(OneMHR(a = a, b = b, c = c+1e-09, d = d)))))  
    }
  }
  else { # all regular cases with reasonable values across the confusion matrix
    return (list(
                 (logFA(a = a, b = b, c = c, d = d) - logHR(a = a, b = b, c = c, d = d) - log(OneMFA(a = a, b = b, c = c, d = d)) + log(OneMHR(a = a, b = b, c = c, d = d)))
                 /
                 (logFA(a = a, b = b, c = c, d = d) + logHR(a = a, b = b, c = c, d = d) + log(OneMFA(a = a, b = b, c = c, d = d)) + log(OneMHR(a = a, b = b, c = c, d = d)))
            ))  
  }
}


# Kappa
kappa=function(a,b,c,d){
  n=a+b+c+d
  out=((a+d)/n-((a+b)*(a+c)+(c+d)*(d+b))/n^2)/(1-((a+b)*(a+c)+(c+d)*(d+b))/n^2)
  return(out)
}
wppp=function(a,b,c,d){
  a^2/((a+b) * (a+c))
}

ed=function(a,b,c,d){
  
  # Definitions
  N=a+b+c+d
  p=(a+c)/N
  pa=a/(a+b)
  pb=c/(c+d)
  Tres=a+b
  
  Dnull=-2*N*(p*log(p)+(1-p)*log(1-p))
  # print(Dnull)
  Dresid=-2*(Tres*(pa*log(pa)+(1-pa)*log(1-pa))+((N-Tres)*(pb*log(pb)+(1-pb)*log(1-pb))))
  # print(Dresid)
  Dexpl=(Dnull-Dresid)/Dnull
  return(Dexpl)
  
}

# combine them in a function
nsdm.all.metrics=function(a,b,c,d){
  out=c(sens(a,c),
        spec(b,d),
        acc(a,b,c,d),
        ppv(a,b),
        npv(c,d),
        jaccard(a,b,c),
        tss(a,b,c,d),
        kappa(a,b,c,d),
		sedi(a,b,c,d))
  names(out)=c("Sensitivity","Specificity","Accuracy","PPV","NPV","Jaccard","TSS","Kappa","SEDI")
  return(out)
}
