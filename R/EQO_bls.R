#' Ensemble Quotient Optimization with a boolean least square (BLS) algorithm
#' 
#' This function performs BLS-based Ensemble Quotient Optimization 
#' towards annotation-free microbiome coarse-graining.  
#' 
#' @param M A matrix or dataframe of taxa in the microbiome (samples as rows and taxa as columns).
#' @param y A vector of a continuous functional trait (length as number of samples). 
#' @param Nmax A number (default 10) specifying the maximal number of taxa allowed in the final group.
#' @param K A sufficiently large number (default 100) required for linearization by the optimizer. A smaller value of K would favor faster optimization by forcing the relative abundance of the selected group to be higher.
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
#' To run BLS more effciently, it is strongly recommended to reduce the total number of taxa
#' by removing the rare singletons (e.g., < 200 taxa in the final input) and apply  
#' strong regularization (e.g., Nmax<10).
#' 
#' @seealso 
#' [EQO_ga],
#' [CAN],
#' 
#' @references 
#' [https://www.biorxiv.org/content/10.1101/2022.08.02.502537v1]
#' 
#' @export
#' 
#' @examples 
#' # This may take several minutes.
#' EQO_bls(Microbiome,trait)

EQO_bls<-function(M,y,Nmax=5,K=100){
	
	# if(missing(Nmax)){Nmax<-5}
	# if(missing(K)){K<-100}
	
  M = as.matrix(M)
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
	
	result <- functionWithgurobi(model)
	
	members<-which(result$x[-1][1:n]>0)
		
	x<-which(result$x[-1][1:n]>0)
	members<-colnames(M)[x]
	abundance<-rowSums(cbind(rep(0,m),rep(0,m),M[,x]))
	R2<-cor(abundance,yn)^2
	return(list(x=x,members=members,abundance=abundance,y=R2))
}

functionWithgurobi<-function(model) {
   if (!requireNamespace("gurobi", quietly = TRUE)) {
        warning("The gurobi package must be installed to use this functionality")
        return(NULL)
    }
    gurobi::gurobi(model)
}



