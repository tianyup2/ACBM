rm(list=ls())
library(Rcpp)
library(ACBM)

log_Vnt <- function(n,m,gamma_mfm,lambda_mfm){
  # assume gamma=1, t ranges from 1 to n+1
  vnt_list=rep(0,n+m)
  nTrunc=1e2
  for(t in 1:(n+m)){
    r=-Inf
    for(k in t:(t+nTrunc-1)){
      b=lgamma(k+1)-lgamma(k-t+1)-(lgamma(n+k*gamma_mfm)-lgamma(k*gamma_mfm))+dpois(k-1,lambda_mfm,log = TRUE)
      m=max(b,r)
      r=log(exp(r-m)+exp(b-m))+m
    }
    vnt_list[t]=r
  }
  return(vnt_list)
}

log_Vnt_Cons <- function(n,K,gamma_mfm,lambda_mfm){
  # K is the upper bound of components
  dpList <- dpois(1:K,lambda = lambda_mfm)
  dpList <- log(dpList/sum(dpList))
  
  vnt_list=rep(0,K)
  for(t in 1:K){
    r <- -Inf
    for(k in t:K){
      b <- lgamma(k+1)-lgamma(k-t+1) - (lgamma(n+k*gamma_mfm)-lgamma(k*gamma_mfm)) + dpList[k]
      m <- max(b,r)
      r <- log(exp(r-m)+exp(b-m))+m
    }
    vnt_list[t] <- r
  }
  return(vnt_list)
}

getMarginal <- function(y_mat, log_Vnt_list, gamma_mfm, n_ml, a0, b0){
  # if marginal likelihood is correct:
  col_ml <- c()
  m=dim(y_mat)[2]
  pb <- txtProgressBar(min = 0, max = m, style = 3)
  for(i in 1:m){
    col_ml[i]<-marginal_column_bin(y = y_mat[,i], a0 = a0, b0 = b0, 
                                   log_Vnt = log_Vnt_list, gamma_mfm = gamma_mfm, n_rep = n_ml)
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(col_ml)
}

getBiClst <- function(y_mat, n, m, col_ml, 
                      gamma_mfm, log_Vnt, log_Vnt_Cons_K, col_assign, row_assign,
                      col_size, row_size, col_n_clst, row_n_clst, n_iter,
                      n_rep=1e2, a0, b0, ifProg=TRUE){
  
  History <- vector("list", n_iter)
  if(ifProg){
    pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  }
  for(i_iter in 1:n_iter)
  {
    consBinClstK(y_mat = y_mat, col_ml = col_ml, 
                 log_Vnt = log_Vnt, log_Vnt_Cons_K = log_Vnt_Cons_K,
                 col_assign = col_assign, row_assign =  row_assign,
                 col_size = col_size, row_size =  row_size,
                 col_n_clst =  col_n_clst, row_n_clst =  row_n_clst,
                 a0 =  a0, b0 =  b0,
                 n_iter = 1, n_rep =  n_rep, gamma_mfm =  gamma_mfm)
    
    
    col_assign_temp <- duplicate(col_assign, shallow = TRUE)
    col_size_temp <- duplicate(col_size, shallow = TRUE)
    col_n_clst_temp <- duplicate(col_n_clst, shallow = TRUE)
    
    row_assign_temp <- duplicate(row_assign, shallow = TRUE)
    row_size_temp <- duplicate(row_size, shallow = TRUE)
    row_n_clst_temp <- duplicate(row_n_clst, shallow = TRUE)
    
    # initialize row cluster
    
    
    History[[i_iter]] <- list(col_assign=col_assign_temp,row_assign=row_assign_temp,
                              col_size=col_size_temp,row_size=row_size_temp,
                              col_n_clst=col_n_clst_temp, row_n_clst=row_n_clst_temp)
    
    if(ifProg){
      setTxtProgressBar(pb, i_iter)
    }
  }
  if(ifProg){
    close(pb)
  }
  return(History)
}


getDahl_colClst <- function(Fit, burn)
{
  iters <- Fit[-(1:burn)]
  niters <- length(iters)
  membershipMatrices <- lapply(iters, function(x){
    col_assign <- x$col_assign
    outer(col_assign, col_assign, FUN = "==")
  })
  membershipAverage <- Reduce("+", membershipMatrices)/niters
  SqError <- sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                    av = membershipAverage)
  DahlIndex <- which.min(SqError)
  DahlAns <- iters[[DahlIndex]]
  attr(DahlAns, "iterIndex") <- burn + DahlIndex
  attr(DahlAns, "burnin") <- burn
  return(DahlAns)
}

getDahl_rowClst <- function(Fit, burn, col_Dahl){
  iters <- Fit[-(1:burn)]
  niters <- length(iters)
  len_assign <- length(iters[[1]]$row_assign[,1])
  m <- length(col_Dahl)
  # initialize membership matrix
  
  membershipMatrix <- list()
  
  for(j in 1:m){
    membershipMatrix[[j]] <- matrix(0,nrow = len_assign,ncol = len_assign)
  }
  
  # calculate the expected correlation matrix
  count <- 0
  for(i in 1:niters){
    if(fossil::rand.index(iters[[i]]$col_assign, col_Dahl)==1){
      for(j in 1:m){
        row_assign_temp <- iters[[i]]$row_assign[,j]
        membershipMatrix[[j]] <- membershipMatrix[[j]] + outer(row_assign_temp,row_assign_temp, FUN = "==")
      }
      count <- count + 1
    }
  }
  
  for(j in 1:m){
    membershipMatrix[[j]] <- membershipMatrix[[j]]/count    
  }
  
  SqError <- c()
  accMax <- 0
  for(i in 1:niters){
    if(fossil::rand.index(iters[[i]]$col_assign, col_Dahl)==1){
      accLoss <- 0
      for(j in 1:m){
        row_assign_temp <- iters[[i]]$row_assign[,col_Dahl[j]]
        membershipMatrix_temp <- outer(row_assign_temp,row_assign_temp, FUN = "==")
        accLoss <- accLoss + sum((membershipMatrix_temp-membershipMatrix[[j]])^2)
      }
      if(accLoss > accMax){
        accMax <- accLoss
      }
      SqError[i] <- accLoss
    }
  }
  
  SqError[is.na(SqError)] <- accMax
  
  DahlIndex <- which.min(SqError)
  DahlAns <- iters[[DahlIndex]]
  attr(DahlAns, "iterIndex") <- burn + DahlIndex
  attr(DahlAns, "burnin") <- burn
  return(DahlAns)
}

genData <- function(schm, n, m, seed = 1){
  set.seed(seed)
  y_mat <- matrix(0,nrow = n,ncol = m)
  rowList <- list()
  if(schm == 1){
    # assume if schema = 1, the model is generated under our model
    # assume m = 20 and 80, n = 100, 300, 1000
    # assume three huge cluster and three trivial cluster
    # each trivial cluster has only one question, and only one student cluster
    # if m = 20, question cluster size: 3,2,2,1,1,1
    # if m = 60, question cluster size: 5,3,3,1,1
    if(m == 20){
      colIdx <- c(rep(1,7),rep(2,6),rep(3,5),4,5)
      colSize <- c(3,2,2,1,1)
      param <- list()
      param[[1]] <- seq(0.1,0.9,length.out = 3)
      param[[2]] <- seq(0.1,0.6,length.out = 2)
      param[[3]] <- seq(0.4,0.9,length.out = 2)
      param[[4]] <- c(0.1)
      param[[5]] <- c(0.9)
      
      pMat <- matrix(0,n,m)
      
      for(j in 1:5){
        idxTemp <- which(colIdx == j)
        idxLen <- length(idxTemp)
        idxRow <- c(1:colSize[j],sample(colSize[j],n-colSize[j],replace = TRUE,prob = rep(1,colSize[j])))
        for(i in 1:n){
          y_mat[i,idxTemp] <- rbinom(idxLen,1,prob = param[[j]][idxRow[i]])
          pMat[i,idxTemp] <- param[[j]][idxRow[i]]
        }
        rowList[[j]] <- idxRow
      }
      
    }else if(m == 60){
      colIdx <- c(rep(1,21),rep(2,20),rep(3,17),4,5)
      colSize <- c(4,3,3,1,1)
      param <- list()
      param[[1]] <- seq(0.1,0.9,length.out = 4)
      param[[2]] <- seq(0.1,0.6,length.out = 3)
      param[[3]] <- seq(0.4,0.9,length.out = 3)
      param[[4]] <- c(0.1)
      param[[5]] <- c(0.9)
      
      pMat <- matrix(0,n,m)
      
      for(j in 1:5){
        idxTemp <- which(colIdx == j)
        idxLen <- length(idxTemp)
        idxRow <- c(1:colSize[j],sample(colSize[j],n-colSize[j],replace = TRUE,prob = rep(1,colSize[j])))
        for(i in 1:n){
          y_mat[i,idxTemp] <- rbinom(idxLen,1,prob = param[[j]][idxRow[i]])
          pMat[i,idxTemp] <- param[[j]][idxRow[i]]
        }
        rowList[[j]] <- idxRow
      }
      
    }else{
      print("invalid column number")
    }
  }else{
    # assume if schema not equal to 1, the model is generated under Rasch model
    # assume m = 20 and 80, n = 100, 300, 1000
    # assume three huge cluster and three trivial cluster
    # each trivial cluster has only one question, and only one student cluster
    # if m = 20, question cluster size: 6,6,8
    # if m = 80, question cluster size: 24, 26, 30
    # colParam: -2, 0, 2
    # rowParam: -2, 0, 2
    if(m == 20){
      colIdx <- c(rep(1,10),rep(2,10))
      rowLabel <- c(1:3,sample(1:3,n-3,replace = TRUE,prob = rep(1,3)))
      paramC <- c(-0.5,0.5)
      paramR <- c(-2,0,2)
      paramCol <- rep(0,m)
      paramRow <- rep(0,n)
      for(j in 1:2){
        idxCol <- which(colIdx == j)
        idxColLen <- length(idxCol)
        # paramCol[idxCol] <- rnorm(idxColLen,mean = param[j],sd = 0.1)
        paramCol[idxCol] <- paramC[j]
      }
      
      for(j in 1:3){
        idxRow <- which(rowLabel == j)
        idxRowLen <- length(idxRow)
        # paramRow[idxRow] <- rnorm(idxRowLen,mean = param[j],sd = 0.1)
        paramRow[idxRow] <- paramR[j]
        
        rowList[[j]] <- rowLabel
      }
      
      pMat <- matrix(0,n,m)
      for(j in 1:m){
        for(i in 1:n){
          probTemp <- exp(paramRow[i]-paramCol[j])/(1+exp(paramRow[i]-paramCol[j]))
          pMat[i,j] <- probTemp
          y_mat[i,j] <- rbinom(1,1,prob = probTemp)
        }
      }
      
    }else if(m == 60){
      colIdx <- c(rep(1,30),rep(2,30))
      rowLabel <- c(1:3,sample(1:3,n-3,replace = TRUE,prob = rep(1,3)))
      
      paramC <- c(-0.5,0.5)
      paramR <- c(-2,0,2)
      
      paramCol <- rep(0,m)
      paramRow <- rep(0,n)
      for(j in 1:2){
        idxCol <- which(colIdx == j)
        idxColLen <- length(idxCol)
        # paramCol[idxCol] <- rnorm(idxColLen,mean = param[j],sd = 0.1)
        paramCol[idxCol] <- paramC[j]
      }
      
      for(j in 1:3){
        idxRow <- which(rowLabel == j)
        idxRowLen <- length(idxRow)
        # paramRow[idxRow] <- rnorm(idxRowLen,mean = param[j],sd = 0.1)
        paramRow[idxRow] <- paramR[j]
        
        rowList[[j]] <- rowLabel
      }
      
      pMat <- matrix(0,n,m)
      for(j in 1:m){
        for(i in 1:n){
          probTemp <- exp(paramRow[i]-paramCol[j])/(1+exp(paramRow[i]-paramCol[j]))
          pMat[i,j] <- probTemp
          y_mat[i,j] <- rbinom(1,1,prob = probTemp)
        }
      }
      
    }else{
      print("invalid column number")
    }
  }
  out <- list()
  out$data <- y_mat
  out$colIdx <- colIdx
  out$rowIdx <- rowList
  out$pmat <- pMat
  return(out)
}
n <- 8e2
m <- 20

yout <- genData(schm = 1,n,m,seed = 1)

y_mat <- yout$data

gamma_mfm <- 1
lambda_mfm <- 1

log_Vnt_list <- list()

for(i in 1:n){
  log_Vnt_list[[i]] <- log_Vnt(i,0,gamma_mfm = gamma_mfm,lambda_mfm = lambda_mfm)
}

log_Vnt_Cons_K <- list()

mTotal <- (m+1)/2

for(i in 1:as.integer((m+1)/2)){
  log_Vnt_Cons_K[[i]] <- log_Vnt_Cons(n,i,gamma_mfm = gamma_mfm,lambda_mfm = lambda_mfm)
}

set.seed(1)
a0 <- 1e-2
b0 <- 1e-2
# initC <- length(unique(yout$colIdx))
initC <- 3
n_rep <- 2e2
n_iter <- 2e2
n_ml <- 1e4

col_assign <- rep(0,m)
col_size <- rep(0,m)

# col_assign[1:m] <- c(yout$colIdx)
col_assign[1:m] <- as.integer(c(1:initC,sample(1:(initC), size = m - initC, replace = TRUE)))
col_size[1:initC] <- table(col_assign)
col_n_clst <- c(initC)

rAssign <- list()
rSize <- list()
rN <- c()

row_assign <- matrix(0,nrow = n,ncol = m)
row_size <- matrix(0,nrow = n,ncol = m)
row_n_clst <- rep(0,m)

for(i in 1:col_n_clst){
  
  initC <- as.integer((col_size[i]+1)/2)

  rAssign[[i]] <- as.integer(c(1:initC,sample(1:(initC), size = n - initC, replace = TRUE)))
  # rAssign[[i]] <- clstLabelList[[i]]
  
  # rAssign[[i]] <- c(yout$rowIdx[[i]])
  
  row_assign[,i] <- rAssign[[i]]
  
  rSize[[i]] <- as.integer(table(rAssign[[i]]))
  
  rN[i] <- length(rSize[[i]])
  row_n_clst[i] <- rN[i]
  
  row_size[1:row_n_clst[i],i] <- as.integer(rSize[[i]])
}

# col_ml <- getMarginal(y_mat = y_mat, log_Vnt = log_Vnt_list, 
#                       gamma_mfm = gamma_mfm, n_ml = n_ml,
#                       a0 = a0, b0 = b0)

col_ml <- rep(0,m)

for(i in 1:m){
  y_sum <- sum(y_mat[,i])
  
  atemp <- a0 + y_sum
  btemp <- b0 + n - y_sum
  col_ml[i] <- (lgamma(atemp) + lgamma(btemp) - lgamma(atemp + btemp))-
    (lgamma(a0) + lgamma(b0) - lgamma(a0 + b0))
}

t1 <- Sys.time()

bicluster_result <- getBiClst(y_mat = y_mat,n = n, m = m,col_ml = col_ml,
                              gamma_mfm = gamma_mfm, log_Vnt = log_Vnt_list, log_Vnt_Cons_K = log_Vnt_Cons_K,
                              col_assign = col_assign,row_assign = row_assign,
                              col_size = col_size, row_size = row_size,
                              col_n_clst = col_n_clst,row_n_clst = row_n_clst,n_iter = n_iter,
                              n_rep=n_rep, a0 = a0, b0 = b0)

t2 <- Sys.time()

result <- list()

bicluster_col_Dahl <- getDahl_colClst(bicluster_result,burn = n_iter/2)

result$col_assign <- duplicate(bicluster_col_Dahl$col_assign, shallow = TRUE)

bicluster_row_dahl <- getDahl_rowClst(bicluster_result,burn = n_iter/2,col_Dahl = result$col_assign)

result$row_assign <- duplicate(bicluster_row_dahl$row_assign, shallow = TRUE)
