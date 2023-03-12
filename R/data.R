#' @docType data
#' @name Microbiome
#'
#' @title Demo data to display mEQO package
#' @description A small sample dataset of a biological experiment to display all functions in package "mEQO".  
#'
#' @format Three arrays included.
#' \itemize{
#'   \item Microbiome - A matrix with 20 rows(samples) and 8 columns(species). 
#'   \item trait_d - A 20*2 matrix indicates the category of each sample(20 samples in 2 categories). 
#'   The target \code{y} in \code{patten} "d".
#'   \item trait - A Vector length 20. The target \code{y} in \code{patten} "c".
#' }
#' 
#' @examples 
#' View(Microbiome)
#' can<-CAN("ga_c",Microbiome,trait,maxIter=100)
#' CAN_plot(can)

#' @source [https://github.com/Xiaoyu2425/mEQO]
NULL