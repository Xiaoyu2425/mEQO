#' Ensemble Quotient Optimization with genetic algorithm (GA)
#' 
#' This function performs GA-based Ensemble Quotient Optimization 
#' towards annotation-free microbiome coarse-graining.  
#' 
#' @param pattern A character speficifying the type of trait, including \code{"u"} for a uniform trait, \code{"c"} for a continuous trait and \code{"d"} for a discrete trait.
#' @param M A matrix or dataframe of taxa in the microbiome (samples as rows and taxa as columns). Rownames and colnames should be specified.
#' @param y A vector or matrix of trait. y is not needed when pattern is \code{"u"}. y is a vector when pattern is \code{"c"}. y is a binary matrix when pattern is \code{"d"} whose rows are samples and columns are categories. Samples belonging to a category is assigned 1, otherwise 0. An example for the discrete case can be found in \code{trait_d}.
#' @param pk A binary vector indicating partially known functional group based on a priori knowledge, with 1 for species forced to be included in the targeted group and 0 for the other unknown species. Length of this vector should be equal to the total number of species in the microbiome. EQO will then search on the basis of the provided partial known group without removing the designated species.
#' @param Nmax A number (default the total number of species) for regularization, specifying the maximal number of taxa allowed in the final group.
#' @param amin A number (default 0) specifying the lower bound of average relative abundance of the final group across all samples (only applied to a uniform trait to avoid trivial solutions).
#' @param amax A number (default 1) specifying the upper bound of average relative abundance of the final group across all samples (only applied to a uniform trait to avoid trivial solutions).
#' @param maxIter A number (default 500) for GA package, specifying the number of generations for GA.
#' @param popSize A number (default 200) for GA package, specifying the initial population size for the GA search. 
#' @param parallel Optional argument for GA package (default TRUE), indicating whether parallel computing is needed. See GA documentation for full details.
#' @param monitor Optional argument for GA package (default plot). \code{monitor=plot} indicates that algorithm implementation is visualized as a plot. \code{monitor=gaMonitor} indicates that algorithm implementation is printed as text. \code{monitor=FALSE} indicates the output should be suppressed. 
#' @return A list with the following components.
#' \itemize{
#'   \item fitness - optimized value for the objective function in GA
#'   \item x - A binary vector specifying presence(1)/absence(0) of taxa in the final group
#'   \item members - A character vector for names of taxa selected in the final group
#'   \item abundance - A numeric vector with relative abundance of the final group in each sample.
#'   \item performance - A number indicating the statistical performance of the group. It is equal to the coefficient of variation when the trait is uniform, for which a smaller value implies higher stability. It is equal to Pearson's correlation coefficient when the trait is continuous, for which a larger absolute value implies stronger association. It is equal to the generalized coefficient of determination when the trait is discrete, for which a larger value implies better discrimination across different categories.
#' }
#' 
#' @seealso 
#' [EQO_bls],
#' [CAN],
#' 
#' @export
#' 
#' @references 
#' [https://www.biorxiv.org/content/10.1101/2022.08.02.502537v1]
#' 
#' @examples 
#' # This may take several minutes.
#' EQO_ga("c",Microbiome,trait,maxIter = 100)
#' EQO_ga("d",Microbiome,trait_d,maxIter = 100)
#' EQO_ga("u",Microbiome,maxIter = 100)

EQO_ga<-function(pattern=NULL,M,y=NULL,pk=NULL,Nmax=NULL,amin=0,amax=1,maxIter=500,popSize=200,parallel=TRUE,monitor=plot){
	
	if(missing(y)){y<-1; pattern<-"u"}
	if(missing(pk)){pk<-rep(0,ncol(M))}
	if(missing(Nmax)){Nmax<-ncol(M)}

  M = as.matrix(M)
	an<-colMeans(M)
	m<-nrow(M)
	n<-ncol(M)
	pen<-sqrt(.Machine$double.xmax)
	
	c0<-function(x){ifelse(min(x-pk)<0,1,0)} # partially known group
	c1<-function(x){(rep(1,n) %*% x)-Nmax} # maximal number of taxa
	c2<-function(x){amin-(an %*% x)} # minimal average relative abundance
	c3<-function(x){(an %*% x)-amax} # maximal average relative abundance

	if(pattern=="c"){
		M0<-apply(M,2,function(x){x-mean(x)})
		y0<-y-mean(y)
		P<-t(M0) %*% M0
		Q<-t(M0) %*% y0 %*% t(y0) %*% M0
		Q2<-t(M0) %*% y0
		
		fitness<-function(x){
			(t(x) %*% Q2)/sqrt((t(x) %*% P) %*% x)-(pen*max(c1(x),0))-(pen*c0(x))
		}
	}

	if(pattern=="d"){
		M0<-apply(M,2,function(x){x-mean(x)})
		y0<-y
		L<-matrix(0,nrow=ncol(y0),ncol=ncol(y0))
		diag(L)<-1/sqrt(colSums(y0))
		P<-t(M0) %*% M0
		Q<-t(M0) %*% y0 %*% L %*% L %*% t(y0) %*% M0
				
		fitness<-function(x){
			((t(x) %*% Q) %*% x)/((t(x) %*% P) %*% x)-(pen*max(c1(x),0))-(pen*c0(x))
		}
	}

	if(pattern=="u"){
		M0<-M
		e<-as.matrix(rep(1,m),ncol=1)
		P<-t(M0) %*% M0 - ((2/n) * (t(M0) %*% e %*% t(e) %*% M0)) + ((1/(n^2)) * (t(M0) %*% e %*% t(e) %*% e %*% t(e) %*% M0))
		Q<-((1/(n^2)) * t(M0) %*% e %*% t(e) %*% M0)
		
		fitness<-function(x){
			((t(x) %*% Q) %*% x)/((t(x) %*% P) %*% x)-(pen*max(c1(x),0))-(pen*max(c2(x),0))-(pen*max(c3(x),0))-(pen*c0(x))
		}
		
	}

	GA<-GA::ga("binary",
		fitness=fitness,
		nBits=n,
		maxiter=maxIter,
		popSize=popSize,
		names=colnames(M0),
		monitor=monitor,
		parallel=parallel)
	
	if(pattern=="c"){
		x<-as.numeric(GA@solution)
		fitness<-max(GA@fitness[!is.na(GA@fitness)])
		members<-colnames(M0)[GA@solution==1]
		abundance<-rowSums(cbind(rep(0,m),rep(0,m),M[,GA@solution==1]))
		r<-cor(abundance,y0)
		return(list(fitness=fitness,x=x,members=members,abundance=abundance,performance=r))
	}

	if(pattern=="d"){
	  solution1 = GA@solution[1,]
	  x<-as.numeric(solution1)
	  fitness<-max(GA@fitness[!is.na(GA@fitness)])
	  members<-colnames(M)[solution1==1]
	  abundance<-rowSums(cbind(rep(0,m),rep(0,m),M[,solution1==1]))
	  s<-abundance-mean(abundance)
	  R2<-(t(s) %*% y0 %*% L %*% L %*% t(y0) %*% s)/(t(s) %*% s)
		return(list(fitness=fitness,x=x,members=members,abundance=abundance,performance=R2))
	}

	if(pattern=="u"){
		x<-as.numeric(GA@solution)
		fitness<-max(GA@fitness[!is.na(GA@fitness)])
		members<-colnames(M0)[GA@solution==1]
		abundance<-rowSums(cbind(rep(0,m),rep(0,m),M0[,GA@solution==1]))
		CV<-sd(abundance)/mean(abundance)
		return(list(fitness=fitness,x=x,members=members,abundance=abundance,performance=CV))
	}
}




