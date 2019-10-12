# This file will contain function-ized versions of the ANALYTICS code. (add Gibbs?)

# These can then be zipped up into a package and installed/used directly from GitHub with:
# devtools::install_github("IsmaelSSE/Bayesian_ReplicationProject")

# (though we probably won't get as good as this guy: https://github.com/kthohr/BMR)
library(magrittr)
library(dplyr)
library(ggplot2)


#This function replicates Korobolis' "mlag2" function:
#custom function takes matrix x and lags it y times
mlag2 <- function(x,y){
  x %<>% as.matrix()
  #create row of 0s, length of matrix to lag
  emptyrow <- rep(0, ncol(x)) %>% as.data.frame() %>% t()
  #list of matrices to store all the lags
  mats <- rep(list(x),y)
  #lags of different lengths for each matrix in list
  for (i in 1:y){
    #inserts one row on top
    zeros <- do.call("rbind", replicate(i, emptyrow, simplify = FALSE))
    mats[[i]] <- rbind(zeros, x)
    mats[[i]] <- mats[[i]][1:nrow(x),] #keeps original length
  }
  outt <- Reduce(cbind, mats)
  return(outt)
}


#replicate "wish.m" from Korobolis:
wish <- function(h,n){
  # h = m x m scale matrix
  # n = scalar degrees freedom
  A <- t(chol(h)) %*% matrix(rnorm(nrow(h)*n),nrow=nrow(h),ncol=n)
  return(A %*% t(A))
} #WORKS!


#Interactive now:


#Interactive!:
bvarAnalytics <- function(){

  #Ask for variables:
  path <- readline(prompt = "File path of data and filename:")
  constant <- readline(prompt = "Include intercept? (TRUE/FALSE)") %>% as.logical()
  p <- readline(prompt = "What are the number of lags?") %>% as.numeric()
  prior <- readline(prompt = "Select a prior (1=noninformative, 2=Minnesota, 3=Nat Conj") %>% as.numeric()
  h <- readline(prompt = "How many periods ahead to predict?") %>% as.numeric()
  forecast_method <- readline(prompt = "Direct forecasts (0) or iterated forecasts (1)?") %>% as.numeric()

  #Loading data:
  Yraw <- read.table(path, header=F)
  #Yraw <- read.table("https://raw.githubusercontent.com/IsmaelSSE/Bayesian_ReplicationProject/master/Yraw.dat?token=AGCJYFYSXUYLIEC7TOHHYTC5TN2HA", header=F)

  #Add step to make sure inputs are valid? (if we care..)

  ##create preliminaries
  #True for intercepts, FALSE else
  #constant <- TRUE
  # number of lags
  #p <- 3
  # predict h steps ahead (TRUE), or don't (FALSE)
  forecasting <- ifelse(h==0, FALSE, TRUE)
  # 0 for direct forecasts, 1 for iterated forecasts
  #forecast_method <- 0
  # how many periods ahead to predict
  #h <- 4

  # set a prior with BVAR model
  #prior <- 3


  #DATA HANDLING
  #initial dimensions
  M <- length(Yraw) #how many variables
  Traw <- nrow(Yraw)


  if (forecasting == 1){
    if (h <= 0){
      print("Wrong forecast horizon set")
    }
    if (forecast_method == 0){
      Y1 <- Yraw[h+1:nrow(Yraw),]
      Y2 <- Yraw[2:(nrow(Yraw)-h),]
      Traw = Traw - h - 1
    } else if (forecast_method == 1){
      Y1 <- Yraw
      Y2 <- Yraw
    } else {
      print("Wrong choice of forecast_method")
    }
  } else {
    Y1 <- Yraw
    Y2 <- Yraw
  }


  #This function replicates Korobolis' "mlag2" function:
  #custom function takes matrix x and lags it y times
  mlag2 <- function(x,y){
    x %<>% as.matrix()
    #create row of 0s, length of matrix to lag
    emptyrow <- rep(0, ncol(x)) %>% as.data.frame() %>% t()
    #list of matrices to store all the lags
    mats <- rep(list(x),y)
    #lags of different lengths for each matrix in list
    for (i in 1:y){
      #inserts one row on top
      zeros <- do.call("rbind", replicate(i, emptyrow, simplify = FALSE))
      mats[[i]] <- rbind(zeros, x)
      mats[[i]] <- mats[[i]][1:nrow(x),] #keeps original length
    }
    outt <- Reduce(cbind, mats)
    return(outt)
  }

  # Create LAG matrix (part of X)
  Ylag <- mlag2(Y2,p)



  #Define matrix X that has all RHS variables:
  if (constant){
    #adds a column of ones
    X1 <- cbind(1, Ylag[(p+1):Traw,])
  } else{
    X1 <- Ylag[(p+1):Traw,]
  }


  #Final dimensions
  K <- ncol(X1) #how many cols
  Traw3 <- nrow(X1)


  #make Block diagonoal matrix Z
  Z1 <- kronecker(diag(M),X1)

  #Make Y matrix
  Y1 <- Y1[(p+1):Traw,]

  #Number of actual time periods "T"
  Tp <- Traw - p

  ## Set up forecast:
  if (forecasting){
    if (forecast_method == 0){
      Y <- Y1[1:(nrow(Y1)-1),]
      X <- X1[1:(nrow(X1)-1),]
      Z <- kronecker(diag(M),X)
      Tp <- Tp - 1
    } else {
      Y <- Y1[1:(nrow(Y1)-h),]
      X <- X1[1:(nrow(X1)-h),]
      Z <- kronecker(diag(M),X)
      Tp <- Tp - h
    }
  }


  ##Priors: OLS
  #OLS estimators:
  Y %<>% as.matrix()
  A_OLS <- solve(t(X) %*% X) %*% (t(X) %*% Y)
  #a_OLS <- A_OLS(:)   #THIS WEIRD STEP to merge columns into one column... why?
  a_OLS <- as.vector(A_OLS)

  SSE <- t(Y - X %*% A_OLS) %*% (Y - X %*% A_OLS)
  SIGMA_OLS <- SSE / (Tp-K)

  ##Priors: Hyperparameters

  if (prior == 1){
    #Uninformative Prior -- do nothing
  } else if (prior == 2){
    #Minnesota Prior
    A_prior <- matrix(0, nrow=K, ncol=M)
    a_prior <- as.vector(A_prior)

    #Define Hyperparameters:
    a_bar_1 <- 0.5 #tightness (0 - prior imposed exactly)
    a_bar_2 <- 0.5 #std dev of prior on lags
    a_bar_3 <- 100 #decay

    #list of numbers to store residual vars
    sigma_sq <- rep(0, M)

    #Simple AR run for each var "M" separately
    for (i in 1:M){
      #Create lags of dep var in i-th equation
      Ylag_i <- matrix(mlag2(Yraw[,i],p))
      Ylag_i <- Ylag_i[(p+1):Traw,]
      #Dependent variable in i-th equation
      Y_i <- as.matrix(Yraw[(p+1):Traw,i])
      #OLS estimates of i-the equation ("solve" = inverse)
      alpha_i <- solve(t(Ylag_i)%*%Ylag_i) %*% (t(Ylag_i)%*%Y_i)
      sigma_sq[i] <- 1/(Tp-p+1) * t(Y_i - Ylag_i%*%alpha_i) %*% (Y_i - Ylag_i %*% alpha_i)
    }

    #V = covariance matrix for the prior
    V_i <- matrix(0, nrow=K, ncol=M)

    #index in each equation where self-lags are:
    ind <- matrix(0, nrow = M, ncol = p)
    for (i in 1:M){
      ind[i,] <- cbind(t(as.matrix(seq(constant + i,K,M))))
    }

    #Calculate (8) on page 6: PRIOR COV MATRIX
    for (i in 1:M){
      for (j in 1:K){
        if (constant){
          if(j==1){
            #part 3
            V_i[j,i] <- a_bar_3 * sigma_sq[i]
          } else if (j %in% ind[i,]){
            #part 1
            V_i[j,i] <- a_bar_1 / p^2  #THIS NOT RIGHT???
          } else {
            #part 2
            for (kj in 1:M){
              if (j %in% ind[kj,]){
                ll <- kj
              }
            }
            V_i[j,i] <- (a_bar_2*sigma_sq[i]) / (p^2 * sigma_sq[ll])
          }
        } else {
          #No constant -- no intercept terms
          if (j %in% ind[i,]){
            V_i[j,i] <- a_bar_1 / p^2
          } else {
            for (kj in 1:M){
              if (j %in% ind[kj,]){
                ll <- kj
              }
            }
            V_i[j,i] <- (a_bar_2*sigma_sq[i]) / (p^2 * sigma_sq[ll])
          }
        }
      }

    }

    #Make diagonal matrix out of vectorized V_i:
    V_prior <- diag(as.vector(V_i))
    SIGMA <- SIGMA_OLS

  } else if (prior == 3){
    ##Normal-Wishart / Natural Conjugate
    #Hyperparameters on a ~ N(a_prior, SIGMA x V_prior)
    A_prior <- matrix(0, nrow=K, ncol=M)
    a_prior <- as.vector(A_prior)
    V_prior <- 10*diag(K)

    #Hyperparameters on inv(SIGMA) ~ W(v_prior, )
    v_prior <- M
    S_prior <- diag(M)
    inv_S_prior <- solve(S_prior)
  }


  ##Posterior hyperparameters
  if (prior == 1){
    #Posterior of alpha|Data ~ Multi-T(kron(SSE,inv(X'X)),alpha_OLS,T-K)
    V_post <- solve(t(X) %*% X)
    a_post <- a_OLS
    A_post <- matrix(a_post, nrow = K, ncol = M, byrow = FALSE) #reshape 21x1 into 7x3

    #posterior of SIGMA|Data ~ inv-Wishart(SSE,T-K)
    S_post <- SSE
    v_post <- Tp -K

    #mean and variance of the Multi-t marginal posterior of alpha
    alpha_mean <- a_post
    alpha_var <- (1/(v_post - M - 1)) * kronecker(S_post,V_post)
  } else if (prior == 2){
    #required quantities for the posterior:
    V_post <- solve( solve(V_prior) + kronecker(solve(SIGMA),t(X) %*% X) )
    a_post <- V_post %*% ( solve(V_prior) %*% a_prior + kronecker(solve(SIGMA),t(X)%*%X) %*% a_OLS)
    A_post <- matrix(a_post, nrow = K, ncol = M, byrow = FALSE)
    #mean is just a_post and variance is V_post
    alpha_mean <- a_post
  } else if (prior == 3){
    #required quantities for the posterior
    #for alpha:
    V_post <- solve( solve(V_prior) + t(X) %*% X )
    A_post <- V_post %*% ( solve(V_prior) %*% A_prior + t(X) %*% X %*% A_OLS )
    a_post <- as.vector(A_post)

    #For SIGMA
    S_post <- SSE + S_prior + t(A_OLS) %*% t(X) %*% X %*% A_OLS + t(A_prior) %*% solve(V_prior) %*% A_prior - t(A_post) %*% ( solve(V_prior) + t(X) %*% X) %*% A_post
    v_post <- Tp + v_prior

    #mean and var of Multi-t marginal
    alpha_mean <- a_post
    alpha_var <- (1/(v_post - M - 1)) * kronecker(S_post,V_post)
  }



  ##Predictive Inference
  #Note: Korobolis' matlab code doesn't allow intercept (does this?)
  X_tplus1 <- if (constant) {if(p < 2){
    t(matrix(c(1,Y[Tp,])))
  } else {
    t(matrix(c(1,Y[Tp,],X[Tp,2:(M*(p-1)+1)])))
  }
  } else {
    if(p < 2){
      t(matrix(c(Y[Tp,])))
    } else {
      t(matrix(c(Y[Tp,],X[Tp,2:(M*(p-1)+1)])))
    }
  }
  #Hacky fix to prevent backwards-indexing, e.g. when p=1, then X[,(2:1)] indexes backwards, unlike in MATLAB

  Pred_mean <- X_tplus1 %*% A_post

  ##HOW TO GET PREDICTIONS??

  #put next row (pred mean) into graph:
  Ynext <- rbind(Yraw,Pred_mean)
  colnames(Ynext) <- c("inflation","unemployment","interestRate")
  Ynext$time <-  seq.int(nrow(Ynext))

  # Everything on the same plot
  Ynext <- reshape2::melt(Ynext, "time")
  plotnext <- ggplot(Ynext, aes(x=time, y=value, color=variable)) + geom_line()


  #List of results:
  resultsList <- list("Sigma Squares" = sigma_sq, "Vector" = V_i, "Plot" = plotnext)   #NOTE: sigma_sq only exists for Minnesota
  return(resultsList)
}


#C:/Users/nickp/Desktop/Yraw.dat
#testObj <- bvarAnalytics()
#^ This list BVAR object has however many items you add to "resultsList" above

