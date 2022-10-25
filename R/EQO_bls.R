#' Ensemble Quotient Optimization with a boolean least square (BLS) algorithm
#' 
#' This function performs BLS-based Ensemble Quotient Optimization 
#' towards annotation-free microbiome coarse-graining.  
#' 
#' @param M A matrix of taxa in the microbiome (samples as rows and taxa as columns).
#' @param y A vector of a continuous functional trait (length as number of samples). 
#' @param Nmax A number (default 10) specifying the maximal number of taxa allowed in the final group.
#' @param K A sufficiently large number (default 100) required for linearization by the optimizer.
#' 
#' @return A list with the following components.
#' \itemize{
#'   \item x - A binary numeric vector specifying presence/absence of taxa in the final group
#'   \item members - A character vector Names of taxa selected in the final group.
#'   \item abundance - A numeric vector with relative abundance of the final group in each sample.
#'   \item y - Coefficient of determination for the linear regression. 
#' }
#' 
#' @details The BLS-based algorithm searches for a group of species
#' which as a whole is most strongly correlated with a continuous functional trait variable.    
#' For a uniform or discrete trait variable, or larger-scale problems 
#' (e.g., number of samples > 200 or number of species > 200) please turn to \code{EQO_ga()}.
#' 
#' @export

EQO_bls<-function(M,y,Nmax,K){
	
	if(missing(Nmax)){Nmax<-5}
	if(missing(K)){K<-100}
	
	m<-nrow(M)
	yn<-y/mean(y)
	n<-ncol(M)
	M_aug<-cbind(1,M)
	L= t(M_aug)%*%M_aug
	b= -2*(t(M_aug)%*%yn)
	
	model          <- list()
	model$A        <- rbind(c(rep(0,2*n+1), 1),c(0,rep(0,n),rep(1,n),0))
	model$rhs      <- c(K,Nmax)
	model$sense    <- c('<','<')
	model$obj      <- c(b, rep(0, n+1))
	model$Q        <- matrix(0,nrow=2*n+2, ncol=2*n+2)
	model$Q[1:(n+1),1:(n+1)]<-L
	model$vtype    <- c(rep('C', n+1), rep('B', n), 'C')
	qcs = list()
	for (i in 1:n) {
  		qcs[[i]]       <- list()
  		qcs[[i]]$Qc    <- Matrix::spMatrix(2*n+2, 2*n+2, c(n+i+1), c(2*n+2), c(-1.0))
  		qcs[[i]]$rhs   <- 0
  		qcs[[i]]$q     <- c(rep(0, i), 1, rep(0, 2*n+1-i))
  		qcs[[i]]$sense <- '='
	}
	model$quadcon <- qcs
	
	result <- gurobi::gurobi(model)
	
	members<-which(result$x[-1][1:n]>0)
		
	x<-which(result$x[-1][1:n]>0)
	members<-colnames(M)[x]
	abundance<-rowSums(cbind(rep(0,m),rep(0,m),M[,x]))
	R2<-cor(abundance,yn)^2
	return(list(x=x,members=members,abundance=abundance,y=R2))
}

