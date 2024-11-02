#' @title RE Measure of Explained Variation of the Model
#'
#' @description This package provides tools to evaluate the proportion of variability in survival times explained by covariates using Schempers' RE measure. The measure assesses how well the model explains the variability in the data, which is important for understanding the impact of covariates in survival analysis.
#'
#' The RE measure is sourced from:
#' \url{https://ibmi3.mf.uni-lj.si/ibmi-english/biostat-center/programje/Re.r}
#'
#' The package includes the following functions:
#' 
#' \itemize{
#'   \item \code{re} - Function for calculating the RE measure.
#'   \item \code{summary.re} - Function for calculating the RE measure over time (cumulative and in specified time windows).
#' }
#' 
#' The \code{re} function can be called with:
#' 
#' \code{re(fit)}
#' 
#' It works with various model objects, including \code{coxph}, \code{survreg}, and \code{aareg}.
#'
#' @details
#' The \code{re} function calculates the RE measure, which consists of the following components:
#' - \strong{Re}: The primary RE measure.
#' - \strong{Re.imp}: The RE measure assuming constant coefficients after the last event time.
#' - \strong{Re.fix}: The RE measure assuming a constant value of coefficients throughout the follow-up time. (Note: Returns NA for counting process type times.)
#' - \strong{se}: Standard error of the RE measure.
#' - \strong{sen0}: Standard error calculated with expected Q.
#' - \strong{C}: C-index generalized to allow for time-dependent covariates, ensuring independence from the censoring process.
#' - \strong{r2nw}: An unweighted measure.
#'
#' @note
#' Left-truncated data is not yet implemented. When calculating \code{Re.imp}, the dataset should not be split after the last event time; each data line should represent a different subject.
#'
#' @param fit A survival model object (e.g., \code{coxph}, \code{survreg}, \code{aareg}).
#' @param ... Additional arguments passed to methods.
#'
#' @return A list containing the various RE measures and their associated statistics.
#'
#' @references
#' - Schemper, M. and R. Henderson (2000). Predictive accuracy and explained variation in Cox regression. Biometrics 56, 249--255.
#' @examples
#' # Load necessary libraries
#' library(survival)
#' library(PAmeasures)
#' # Use Mayo Clinic Primary Biliary Cirrhosis Data
#' data(pbc)
#'
#' pbc <- pbc %>% 
#'   filter(!is.na(trt)) %>% 
#'   mutate(log_albumin = log(albumin),
#'          log_bili = log(bili),
#'          log_protime = log(protime),
#'          status = ifelse(status == 2, 1, 0))
#'
#' # Fit a full Cox PH model
#' fit.coxph.full <- coxph(Surv(time, status) ~ age + log_albumin + 
#'                          log_bili + log_protime + edema, 
#'                          data = pbc, x = TRUE, y = TRUE)
#' 
#' # Calculate the RE measure
#' re_result <- pam.re(fit.coxph.full)
#' print(re_result$Re)
#' 
#' # Summarize the RE measure over time
#' summary_result <- pam.summary(re_result)
#' 
#' @export


pam.re <- function (fit, ...) UseMethod("pam.re")


#' @export
pam.re.coxph <- function(fit,Gmat){
  require(survival)
  Y <- fit$y
  f.type <- attr(Y, "type")  
  if (f.type== "right") Y <- cbind(rep(0,nrow(Y)),Y)
  sort.it <- order(Y[,2],-Y[,3]) 
  bx <- fit$linear.predictors[sort.it] - mean(fit$linear.predictors)
  Y <- Y[sort.it,]
  
  n.ind <- nrow(Y)
  ti <- sort(unique(Y[Y[,3]==1,2]))     #ordered unique times of events
  n.times <- length(ti)         #number of unique times of events 
  
  
  r2<-NULL
  if(missing(Gmat)){
    km.inv <- my.survfit(Y[,1],Y[,2],Y[,3])   #now includes time-dependent data
    G <- km.inv$surv.i2[km.inv$n.event != 0]                
    dSt <- 1/G            #is this needed (we already have only at event times)
    Gmat <- matrix(G,byrow=T,ncol=length(G),nrow=n.ind) #times in columns, ind. in rows
  }
  else if(is.null(dim(Gmat))){
    if(length(Gmat)!=n.times)stop("Wrong length of weights")
    else Gmat <- matrix(Gmat,byrow=T,ncol=length(Gmat),nrow=n.ind)
  }
  else if(any(dim(Gmat)!=c(n.ind,n.times)))stop("Wrong dimension of weights")
  Gmat <- Gmat[sort.it,]          #sort with respect to event times
  
  
  ebxun <- ebx <- exp(bx)
  rg <- meanr <- smin <- newvar1 <- newvar2 <- newvar3 <- eden <- enum <-  rep(NA,n.times)
  for(it in 1:length(ti)){
    k <- sum(Y[Y[,3]==1,2]==ti[it])     #number of ties
    
    inx <- (Y[,2] >= ti[it] & Y[,1]< ti[it])  #at risk
    #reverse the order to get the right ranks
    rangi <- (rank(-bx[inx]) - .5)/Gmat[inx,it]+.5  #ranks for each individual
    
    rg[it] <- sum(rangi[1:k]/(Gmat[inx,it])[1:k]) #rank of the one who had the event (sum of ranks if ties)
    
    r0 <- mean(rangi)
    meanr[it] <-  r0*sum(1/((Gmat[inx,it])[1:k])) #mean rank at this time 
    #browser()
    rP <- (.5+ .5*sum((Y[inx,2] == ti[it] & Y[inx,1]< ti[it])/Gmat[inx,it]))
    smin[it] <- rP*sum(1/((Gmat[inx,it])[1:k])) #min rank (1 if no ties)
    
    pi <- ebx[inx]/sum(ebx[inx])      #individual hazard (breslow * ebx)
    newvar1[it] <- sum((rangi-r0)^2*pi/(Gmat[inx,it])^2)  #term1 -  for nom
    newvar2[it] <- sum((rP-r0)^2*pi/(Gmat[inx,it])^2) #term2 -  for denom
    newvar3[it] <- sum((r0-rangi)*(r0-rP)*pi/(Gmat[inx,it])^2)  #term3 - covariance
    enum[it] <- sum((r0-rangi)*pi/(Gmat[inx,it]))   #expected numerator=Q1
    eden[it] <- sum((r0-rP)*pi/(Gmat[inx,it]))  #expected denominator=Q2
    
  }
  num <- sum(meanr-rg)
  den<- sum(meanr-smin)
  eden <- sum(eden)
  enum <- sum(enum)
  r2 <- num/den #the Re measure
  r2nw <- sum((meanr-rg))/sum((meanr-smin))
  newvar <- sqrt(sum(newvar1)/den^2-2*sum(newvar3)*num/den^3+sum(newvar2)*num^2/den^4)
  newvare <- sqrt(sum(newvar1)/eden^2-2*sum(newvar3)*enum/eden^3+sum(newvar2)*enum^2/eden^4)
  
  ebx <- sort(ebx)          #this is for Jfixed
  Jf <- NA
  if(f.type=="right"){
    nd <- length(ebx)
    imat <- matrix(ebx,nrow=nd,ncol=nd)
    jmat <- matrix(ebx,nrow=nd,ncol=nd,byrow=T)
    iimat <- imat/(imat+jmat)
    jjmat <- jmat/(imat+jmat)
    Jf <- 2*sum((jjmat-iimat)[upper.tri(jjmat)])/(sum(upper.tri(jjmat))*2)
  }
  
  #and now for the part after last event time:
  ebxun <- ebxun[Y[,2]>=max(ti)&Y[,3]==0]
  
  km.inv <- my.survfit(Y[,1],Y[,2],Y[,3])   
  nd <- sum(km.inv$n.cens[km.inv$time>max(ti)])
  if((nd>0)&(!any(km.inv$n.cens[km.inv$time>max(ti)]==0))){   #if there is at least 1 indiv. and no time-splitting
    
    ebx <- sort(ebxun)
    sumt <- 0
    for(it in 1:nd){
      sumt1 <- sum((ebx[it]-ebx[1:it])/(ebx[it]+ebx[1:it]))
      sumt <- sumt + sumt1
    }
    Jfin <- sumt/2
    
    Jfin <- (num+Jfin)/(den +  (nd*(nd-1))/4)
  }
  
  else Jfin <- NA     
  
  out <- list(Re=r2,Re.imp=Jfin,Re.fix=Jf,se=newvare,C=(r2+1)/2,sen0=newvar,times=ti,meanr=meanr,ranks=rg,weights=Gmat,perfr=smin,type="coxph",nd=nd,r2nw=r2nw)
  class(out) <- "re"
  out
}

#' @export
pam.re.aareg <- function(fit,Gmat){
  require(survival)
  if(length(fit$y)==0)stop("The aareg model must be fitted with 'y=TRUE' option")
  if(length(fit$x)==0)stop("The aareg model must be fitted with 'x=TRUE' option")
  Y <- fit$y
  sort.it <- order(Y[,2],-Y[,3]) 
  X <- fit$x
  X <- cbind(rep(1,nrow(X)),X)
  Y[Y[,1]==-1,1] <- 0
  Y <- Y[sort.it,]
  X <- X[sort.it,,drop=FALSE]
  n.ind <- nrow(Y)
  
  ti <- sort(unique(Y[Y[,3]==1,2]))     #ordered unique times of events
  n.times <- length(ti)         #number of unique times of events 
  
  #nelson-Aalen:
  d <- survfit(Surv(Y[,1],Y[,2],Y[,3])~1,type="fleming-harrington")
  cumh <- -log(d$surv[d$n.event!=0])
  haz <- c(cumh[1],diff(cumh))
  
  #change the output for ties in surv. time (sum up the contributions in those intervals)
  d1 <- fit$times
  d2 <- c(d1[-1],d1[length(d1)]+1)
  ic <- which(d1!=d2)
  fc <- apply(fit$coef,2,cumsum)
  fc <- fc[ic,]
  fd <- apply(fc,2,function(x)c(x[1],diff(x)))
  #add lines for times at the end 
  
  if(n.times>nrow(fd)){
    nman <- n.times-nrow(fd)
    fd <- rbind(fd,matrix(0,nrow=nman,ncol=ncol(fd)))   
    fd[(n.times-nman+1):n.times,1] <- haz[(n.times-nman+1):n.times]
  }
  
  bx <- (X%*%t(fd))
  
  if(missing(Gmat)){
    km.inv <- my.survfit(Y[,1],Y[,2],Y[,3])   #now includes time-dependent data
    G <- km.inv$surv.i2[km.inv$n.event != 0]              
    dSt <- 1/G
    Gmat <- matrix(G,byrow=T,ncol=length(G),nrow=n.ind) #casi v stolpcih, ljudje v vrsticah 
  }
  else if(is.null(dim(Gmat))){
    if(length(Gmat)!=n.times)stop("Wrong length of weights")
    else Gmat <- matrix(Gmat,byrow=T,ncol=length(Gmat),nrow=n.ind)
  }
  else if(any(dim(Gmat)!=c(n.ind,n.times)))stop("Wrong dimension of weights")
  
  
  rg <- meanr <- smin <- newvar1 <- newvar2 <- newvar3 <- enum <- eden <-  rep(NA,n.times)
  for(it in 1:n.times){
    
    k <- sum(Y[Y[,3]==1,2]==ti[it])     #number of ties
    
    inx <- (Y[,2] >= ti[it] & Y[,1]< ti[it])  #at risk
    pi <- bx[inx,it]      #hazard for each indiv.
    
    rangi <- (rank(-bx[inx,it]) - .5)/Gmat[inx,it]+.5 #ranks for each individual
    
    r0 <- mean(rangi)
    rP <- (.5+ .5*sum((Y[inx,2] == ti[it] & Y[inx,1]< ti[it])/Gmat[inx,it]))
    
    rg[it] <- sum(rangi[1:k]/(Gmat[inx,it])[1:k]) #rank of the one who had the event (sum of ranks if ties)
    meanr[it] <-  r0*sum(1/((Gmat[inx,it])[1:k])) #mean rank at this time 
    smin[it] <- rP*sum(1/((Gmat[inx,it])[1:k]))   
    
    newvar1[it] <- sum((rangi-r0)^2*pi/(Gmat[inx,it])^2)  #term1 -  for nom
    newvar2[it] <- sum((rP-r0)^2*pi/(Gmat[inx,it])^2) #term2 -  for denom
    newvar3[it] <- sum((r0-rangi)*(r0-rP)*pi/(Gmat[inx,it])^2)  #term3 - covariance
    enum[it] <- sum((r0-rangi)*pi/(Gmat[inx,it]))   #expected numerator
    eden[it] <- sum((r0-rP)*pi/(Gmat[inx,it]))  #expected denominator
    
  }
  num <- sum(meanr-rg)
  den<- sum(meanr-smin)
  enum <- sum(enum)
  eden <- sum(eden)
  r2 <- num/den #the Re measure
  newvar <- sqrt(sum(newvar1)/den^2-2*sum(newvar3)*num/den^3+sum(newvar2)*num^2/den^4)
  newvare <- sqrt(sum(newvar1)/eden^2-2*sum(newvar3)*enum/eden^3+sum(newvar2)*enum^2/eden^4)
  out <- list(Re=r2,se=newvare,C=(r2+1)/2,sen0=newvar,type="aareg")
  class(out) <- "re"
  out
}

#' @export
pam.re.survreg <- function(fit,Gmat){
  require(survival)
  Y <- fit$y
  Y <- cbind(rep(0,nrow(Y)),Y)
  if(any(Y[,2]<0))warning("Negative follow-up times")
  #nazaj antilogaritmiramo cas:
  #if(fit$dist%in%c("weibull","exponential","loglogistic","lognormal"))Y[,2] <- exp(Y[,2])
  
  sort.it <- order(Y[,2],-Y[,3]) 
  Y <- Y[sort.it,]
  
  if(fit$dist=="weibull"|fit$dist=="exponential")bx <- -(fit$linear.predictors[sort.it] - fit$coef[1])
  bx <- -fit$linear.predictors[sort.it]
  if(fit$dist=="weibull"|fit$dist=="exponential"){
    gamma <- (1/fit$scale)
    lambda <- exp(-fit$coef[1])^gamma
  }
  
  n.ind <- nrow(Y)
  ti <- sort(unique(Y[Y[,3]==1,2]))     #ordered unique times of events
  n.times <- length(ti)         #number of unique times of events 
  
  r2<-NULL
  
  if(missing(Gmat)){
    km.inv <- my.survfit(Y[,1],Y[,2],Y[,3])   #now includes time-dependent data
    G <- km.inv$surv.i2[km.inv$n.event != 0]              
    dSt <- 1/G
    Gmat <- matrix(G,byrow=T,ncol=length(G),nrow=n.ind) #casi v stolpcih, ljudje v vrsticah 
  } else if(is.null(dim(Gmat))){
    if(length(Gmat)!=n.times)stop("Wrong length of weights")
    else Gmat <- matrix(Gmat,byrow=T,ncol=length(Gmat),nrow=n.ind)
  } else if(any(dim(Gmat)!=c(n.ind,n.times)))stop("Wrong dimension of weights")
  
  
  rg <- meanr <- smin <- newvar1 <- newvar2 <- newvar3 <- eden <- enum  <- rep(NA,n.times)
  for(it in 1:length(ti)){
    k <- sum(Y[Y[,3]==1,2]==ti[it])     #number of ties
    inx <- (Y[,2] >= ti[it] & Y[,1]< ti[it])  #at risk
    if(fit$dist=="weibull"|fit$dist=="exponential"){
      if(it==1) pi <- lambda*exp(bx)[inx]*(ti[it])^gamma            
      else pi <- lambda*exp(bx)[inx]*((ti[it])^gamma-(ti[it-1])^gamma)  
    }
    else if(fit$dist=="gaussian"|fit$dist=="lognormal"){
      bxi <- -bx[inx]   #dam ga nazaj na pozitivne
      casi <- ti
      if(fit$dist=="lognormal") {
        casi <- log(ti)
      }
      S  <- 1-pnorm((casi[it]- bxi)/fit$scale)
      S0 <- 1
      if(it!=1) S0 <- 1-pnorm((casi[it-1]- bxi)/fit$scale)
      pi <- -log(S/S0)
    }     
    else if(fit$dist=="logistic"|fit$dist=="loglogistic"){
      bxi <- -bx[inx]
      casi <- ti
      if(fit$dist=="loglogistic") {
        casi <- log(ti)
      }
      w <- exp((casi[it]- bxi)/fit$scale)
      S <- 1/(1+w)
      S0 <- 1
      if(it!=1){
        w0 <- exp((casi[it-1]- bxi)/fit$scale)
        S0 <- 1/(1+w0)
      }
      pi <- -log(S/S0)
    }     
    rangi <- (rank(-bx[inx]) - .5)/Gmat[inx,it]+.5  #ranks for each individual
    
    r0 <- mean(rangi)
    rP <- (.5+ .5*sum((Y[inx,2] == ti[it] & Y[inx,1]< ti[it])/Gmat[inx,it]))
    
    rg[it] <- sum(rangi[1:k]/(Gmat[inx,it])[1:k]) #rank of the one who had the event (sum of ranks if ties)
    meanr[it] <-  r0*sum(1/((Gmat[inx,it])[1:k])) #mean rank at this time 
    smin[it] <- rP*sum(1/((Gmat[inx,it])[1:k]))   
    
    newvar1[it] <- sum((rangi-r0)^2*pi/(Gmat[inx,it])^2)  #term1 -  for nom
    newvar2[it] <- sum((rP-r0)^2*pi/(Gmat[inx,it])^2) #term2 -  for denom
    newvar3[it] <- sum((r0-rangi)*(r0-rP)*pi/(Gmat[inx,it])^2)  #term3 - covariance
    enum[it] <- sum((r0-rangi)*pi/(Gmat[inx,it]))   #expected numerator
    eden[it] <- sum((r0-rP)*pi/(Gmat[inx,it]))  #expected denominator
    
  }
  num <- sum(meanr-rg)
  den<- sum(meanr-smin)
  enum <- sum(enum)
  eden <- sum(eden)
  r2 <- num/den #the Re measure
  newvar <- sqrt(sum(newvar1)/den^2-2*sum(newvar3)*num/den^3+sum(newvar2)*num^2/den^4)
  newvare <- sqrt(sum(newvar1)/eden^2-2*sum(newvar3)*enum/eden^3+sum(newvar2)*enum^2/eden^4)
  
  out <- list(Re=r2,se=newvare,C=(r2+1)/2,sen0=newvar,times=ti,meanr=meanr,perfr=smin,ranks=rg,weights=Gmat,type="survreg")
  class(out) <- "re"
  out 
}


#' @export
pam.print.re <- function(x, digits=4, ...){
  
  cat("Re measure=",round(x$Re,digits), ",    SE=", round(x$se,digits), "\n \n",sep="")
  if(x$type=="coxph"){
    if(x$nd){
      cat("Number of censored after the last failure time: ", x$nd, "\n")
      if(!is.na(x$Re.imp))cat("Re corrected for the censoring at the end \n (assuming the coeff. and covariate values fixed after the last failure time): ", round(x$Re.imp,digits),"\n \n")
      else cat("Re corrected for the censoring at the end can not be computed - the follow-up \n time after the last event time is split into several intervals \n")
    }
    if(!is.na(x$Re.fix)){
      cat("Integrated version of Re (assuming all coeff. and covariate values fixed):", round(x$Re.fix,digits),"\n")  
    }
  }
  cat("\n")
  cat("Generalized C index=",round(x$C,digits), "\n")
  invisible(x)
}


#' @export
pam.summary.re <- function(object, times, band=5,...){
  
  num <- cumsum(object$meanr-object$ranks)
  den <- cumsum(object$meanr-object$perfr) 
  
  Rti <- num/den
  
  dn <- length(num) - band 
  
  num1 <- num[-(1:band)]
  den1 <- den[-(1:band)]
  
  dRti <- (num1 - num[1:dn])/(den1 - den[1:dn])
  
  dRti <- c(rep(NA,band-1),Rti[band],dRti)
  
  tab <- data.frame(times=object$times,Rti=Rti,dRti=dRti)
  if(!missing(times)) {
    search.time <- function(time,times)max(which(times < time))
    inx <- unlist(lapply(times,search.time, tab$times)  )
    tab[inx,1] <- times
    tab <- tab[inx,]
  }
  tab
}




"my.survfit" <- function(start,stop,event){
  #internal function for calculation of weights (allows counting process notation)
  
  if(missing(start)) start <- rep(0,length(stop))		#if no count. process notation, set everyone to start at 0
  
  data <- data.frame(start=start,stop=stop,event=event)	
  
  data <- data[order(data$start),]			#ordered according to the starting time
  data$Y <- unlist(lapply(data$stop,function(x)sum(data$start<=x)))	#for each time point: the number of indiv. that started before (or at) it
  
  data <- data[order(data$stop,-data$event),]		#order in time, put events before censoring
  
  howmany <- nrow(data)					#total number of rows
  
  lt.t1 <- c(0,data$stop[-nrow(data)])			#aux. vector for searching the stop times	
  lt.t2 <- c(data$stop[-1],data$stop[nrow(data)]+1)
  lt1 <- which(data$stop!=lt.t1)				#first case at a time point
  lt2 <- which(data$stop!=lt.t2)				#last case at a time point
  
  ct.t1 <- c(-1,data$start[-nrow(data)])			#aux. vector for searching the start times
  ct.t2 <- c(data$start[-1],data$start[nrow(data)]+1)
  ct1 <- which(data$start!=ct.t1)				#first case at a time point
  ct2 <- which(data$start!=ct.t2)				#last case at a time point
  n.incom <- rep(0,length(lt1))				#prepare a vector that will contain the no. of incoming indiv. at a time point (length=no. of unique time points)
  inx <- unlist(lapply(data$start[ct1],function(x)min(which(data$stop[lt1]>=x)))) #the index of the first one to fail after a starting point
  incom.sum <- diff(c(0,ct2))				#the number of incoming
  if(length(inx)>1){
    if(inx[1]==inx[2]){					#if indiv. come in at time 0 and first stop time
      incom.sum <- c(sum(incom.sum[1:2]),incom.sum[-(1:2)])    #sum them up
      inx <- unique(inx)					#only one value per each stop time
    }
  }
  n.incom[inx] <- incom.sum				#at each event time, we need the number of patients that got in just before it.
  n.incom2 <- n.incom					
  n.incom2[1] <- n.incom2[1] - sum(data$start==0)		#the number of additional patients to come in at time 0 is set to 0
  
  
  evcum <- cumsum(data$event)			#number of events up to a certain time point
  evcum <- evcum[lt2]				#number of events for t<=t_i
  cencum <- cumsum(data$event==0)			#number of censored
  cencum <- cencum[lt2]				#number of censored t<=t_i
  n.event <- diff(c(0,evcum))			#number of events at each unique event time
  n.cens <- diff(c(-sum(data$start==0),cencum))	#number of censored at each unique event time, the first number is the total of the patients
  n.t.cens <- n.cens - n.incom			#number truly censored: number censored - number of incoming (0 at time 0)
  n.risk <- (data$Y - (0:(howmany-1)))[lt1]-n.incom2 #number at risk at each event time: no.started before or at this time - no. of incoming at this time - 1 per row before this time
  
  kmji <- 1- n.event/n.risk			#proportion of survivors
  
  kmji.i <- 1- n.t.cens/(n.risk-n.event)		#the events are taken out of the risk set before calculating the G(t_i)
  
  km <- (cumprod(kmji))				#km	
  km.i <- cumprod(kmji.i)				#km for censoring			
  
  #correction for the last value (if division by 0) (this happens if the last time is an event time)
  inx <- which(is.na(km.i))				#missing at the end, only needed if not lagged?
  if(length(inx)){
    km.i[inx] <- km.i[min(inx)-1]	# carrie last value forward
  }
  
  km.i2 <- c(1,km.i[-length(km.i)])			#this gives the lagged weight (weight at time just before t)
  
  list(surv=km,surv.i=km.i,n.event=n.event,n.cens=n.t.cens,time=data$stop[lt1],n.risk=n.risk,surv.i2=km.i2)
}
