#' Cross validation-based Aggregation Network (CAN)
#' 
#' This function splits the microbiome dataset into a training set and a test set. 
#' Then it generates a group of taxa by Ensemble Quotient Optimization for the training set  
#' and evaluates the performance of the group with the test set. Finally, the performance of 
#' all taxa selected in any group is aggregated and normalized into its relative importance.
#' 
#' @param method A string speficifying the algorithm to use for EQO. Possible options include \code{"bls_c"} for a continuous variable with BLS, \code{"ga_c"} for a continuous variable with GA, \code{"ga_d"} for a discrete trait with GA and \code{"ga_u"} for a uniform trait with GA.
#' @param M A matrix of taxa in the microbiome (samples as rows and taxa as columns).
#' @param y A vector or matrix of trait. y is optional when pattern is \code{"u"}. y is a vector when pattern is \code{"c"}. y is a binary matrix when pattern is \code{"d"} whose rows are samples and columns are categories.
#' @param fraction A number (default 0.5) indicating the fraction of samples to be randomly selected in the training set.
#' @param tm A number (default 20) indicating the times of cross-validation to perform.
#' @param K (for \code{"bls_c"}) A sufficiently large number (default 100) required for linearization by the optimizer.
#' @param Nmax A number (default 10) for regularization, specifying the maximal number of taxa allowed in the final group.
#' @param pk (for \code{"ga_c"},\code{"ga_d"} or \code{"ga_u"}) A binary vector indicating partially known functional group based on a priori knowledge, with 1 for species forced to be included in the targeted group and 0 for the other unknown species. Length of this vector should be equal to the total number of species in the microbiome. EQO will then search on the basis of the provided partial known group without removing the designated species.
#' @param amin A number (default 0) specifying the lower bound of average relative abundance of the final group across all samples (only applied to a uniform trait to avoid trivial solutions).
#' @param amax A number (default 1) specifying the upper bound of average relative abundance of the final group across all samples (only applied to a uniform trait to avoid trivial solutions).
#' @param maxIter A number (default 500) for GA package, specifying the maximum number of iterations to run before the GA search is halted.
#' @param popSize A number (default 200) for GA package, specifying the initial population size for the GA search. 
#' @param parallel Optional argument for GA package (default TRUE), indicating whether parallel computing is needed. See GA documentation for full details.
#' @param monitor Optional argument for GA package (default FALSE). \code{monitor=plot} indicates that algorithm implementation is visualized as a plot. \code{monitor=gaMonitor} indicates that algorithm implementation is printed as text. \code{monitor=FALSE} indicates the output should be suppressed. 
#' @return A list with the following components.
#' \itemize{
#'   \item output - a list with taxa selected in each time of optimization and the corresponding performances in cross-validation 
#'   \item nodes - a table of nodes for CAN
#'   \item edges - a table of edges for CAN
#' }
#' @export

CAN<-function(method,M,y,fraction,tm,K,pk,Nmax,amin,amax,maxIter=500,popSize=200,parallel=TRUE,monitor=FALSE){
	
	if(missing(fraction)){fraction<-0.5}
	if(missing(tm)){tm<-20}
	if(missing(K)){K<-100}
	if(missing(Nmax)){Nmax<-5}
	if(missing(y)){method<-"ga_u"}
	if(missing(pk)){pk<-rep(0,ncol(M))}
	if(missing(amin)){amin<-0}
	if(missing(amax)){amax<-1}

	size<-round(nrow(M)*fraction,digits=0)
	
	xv.all<-lapply(1:tm,function(t){	
		
		print(paste0(paste0(paste0("Running cross-validation ",t)," out of "),tm))
		
		if(method=="bls_c"){
			seed<-sample(1:nrow(M),size)
			M.train<-M[seed,]
			y.train<-y[seed]
			M.test<-M[-seed,]
			y.test<-y[-seed]
			out<-EQ_bls(M.train,y.train,Nmax,K)
			s<-rowSums(cbind(rep(0,nrow(M.test)),M.test[,which(out$x==1)]))
			y.xv<-cor(s,y.test)
			out.xv<-data.table::data.table(taxa=out$members,perform=y.xv)
			return(out.xv)
		}
		
		if(method=="ga_c"){
			seed<-sample(1:nrow(M),size)
			M.train<-M[seed,]
			y.train<-y[seed]
			M.test<-M[-seed,]
			y.test<-y[-seed]
			out<-EQO_ga("c",M.train,y.train,pk,Nmax,amin,amax,maxIter,popSize,parallel,monitor)
			s<-rowSums(cbind(rep(0,nrow(M.test)),M.test[,which(out$x==1)]))
			y.xv<-cor(s,y.test)
			out.xv<-data.table::data.table(taxa=out$members,perform=y.xv)
			return(out.xv)
		}
		
		if(method=="ga_d"){
			sizes<-round((colSums(y)/sum(colSums(y)))[ncol(y)-1] * size,digits=0)
			sizes<-c(sizes,size-sizes)
			seed<-c()
			for (i in 1:ncol(y)){seed<-c(seed,sample(which(y[,i]==1),sizes[i]))}
			M.train<-M[seed,]
			y.train<-y[seed,]
			M.test<-M[-seed,]
			y.test<-y[-seed,]
			out<-EQO_ga("d",M.train,y.train,pk,Nmax,amin,amax,maxIter,popSize,parallel,monitor)
			s<-rowSums(cbind(rep(0,nrow(M.test)),M.test[,which(out$x==1)]))
			L.test<-matrix(0,nrow=ncol(y.test),ncol=ncol(y.test))
			diag(L.test)<-1/sqrt(colSums(y.test))
			y.xv<-(t(s) %*% y.test %*% L.test %*% L.test %*% t(y.test) %*% s)/(t(s) %*% s)
			out.xv<-data.table::data.table(taxa=out$members,perform=y.xv)
			return(out.xv)
		}
		
		if(method=="ga_u"){
			seed<-sample(1:nrow(M),size)
			M.train<-M[seed,]
			M.test<-M[-seed,]
			y.train<-1
			out<-EQO_ga("u",M.train,y.train,pk,Nmax,amin,amax,maxIter,popSize,parallel,monitor)
			s<-rowSums(cbind(rep(0,nrow(M.test)),M.test[,which(out$x==1)]))
			y.xv<-mean(s)/sd(s)
			out.xv<-data.table::data.table(taxa=out$members,perform=y.xv)
			return(out.xv)
		}
	})
	
	xv.pileup<-lapply(xv.all,function(x){
		pairs<-expand.grid(x$taxa,x$taxa)
		pairs$value<-x$perform
		return(pairs)
	})

	xv.all.dt<-data.table::rbindlist(xv.all)
	nodes<-as.data.frame(xv.all.dt[,sum(perform),by=taxa])
	nodes<-nodes[order(nodes$V1,decreasing=TRUE),]
	colnames(nodes)<-c("taxa","value")
	nodes$value<-nodes$value/max(nodes$value)
	rownames(nodes)<-nodes$taxa
	nodes$id<-rownames(nodes)
	nodes$label<-rownames(nodes)
	
	pileup.dt<-data.table::rbindlist(xv.pileup)
	weight.dt<-plyr::ddply(pileup.dt,~Var1+Var2,dplyr::summarise,weight=sum(value))
	weight.filter<-weight.dt
	weight.mat<-reshape2::dcast(weight.filter,Var1~Var2,value.var="weight")
	rownames(weight.mat)<-as.character(weight.mat$Var1)
	weight.mat<-as.matrix(weight.mat[,-1])
	weight.mat[upper.tri(weight.mat)]<-NA
	diag(weight.mat)<-NA
	edges<-reshape2::melt(weight.mat)
	edges<-edges[!is.na(edges$value),]
	colnames(edges)<-c("from","to","width")
	
	return(list(output=xv.all,nodes=nodes,edges=edges))
}



