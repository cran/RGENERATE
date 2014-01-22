# TODO: Add comment
# 
# Author: ecor
###############################################################################

NULL

#' generate
#' 
#' 
#' 
#' @rdname generate
#' @export

generate <- function (x=NULL,...)  {
	
	
	return(UseMethod("generate",x))
	
}


NULL
#' 
#' 

#' 
#' @rdname generate
#' @method generate default
#' @S3method generate default
#' @aliases generate 
#' @export



generate.default <- function (x,FUN=rnorm,n=100,K=3,names=NULL,cov=NULL,...)  {
	
	if (!is.null(names)) K <- length(names)
	if (!is.null(cov)) K <- ncol(cov)
	
	out <- FUN(n*K,...)
	
	out <- array(out,c(n,K))
	out <- as.data.frame(out)
	if (!is.null(names)) {
		if (length(names)==K) names(out) <- names
	}
# TRASH ....	
#	if (is.null(B)) B <- t(chol(summary(var)$covres))
	
#	if (is.null(xprev)) xprev <- rnorm(ncol(B))
	
	
	if (!is.null(cov)) {
		B <- t(chol(cov))
		out[,] <- ((as.matrix(out)) %*% t(B))
	}
	

	return(out)
	
}


NULL
#' 
#' It generates a multivarite random series according to the model \code{x}
#' 
#' 
#' 
#' 
#' @param x null object or the model used for random generation , e.g. a VAR model as a \code{\link{varest-class}} or \code{\link{varest2-class}} object. Default is \code{NULL}.
#' @param FUN random function of the probability distribution used for noise random generation. Default is \code{\link{rnorm}}. See \url{http://cran.r-project.org/web/views/Distributions.html} 
#' @param n number of generations requested 
#' @param names null object or string vectors or names of the variables to be generated simultaneously. Default is \code{NULL}.
#' @param K number of the variables to be generated simultaneously, i.e. the K parameters of a VAR. It is automatically detected by \code{x}, \code{names} or \code{cov}, if one of these is not \code{NULL}. 
#' @param cov null object or covariance matrix of the random variables to be generated simultaneously. Default is \code{NULL}, not used in case this information can be detected from \code{x}.
#' 
#' @param noise null object or a generic external noise for \code{x} model residuals, e.g. standard white noise, for random generation with the model \code{x}. Default is \code{NULL}. If \code{NULL} the noise is automatically calculated. 
#' @param exogen null object or amatrix or data frame with exogeneous variables (predictors) id requested by \code{x}. Default is \code{NULL}  
#' @param xprev null object or initial condition of the multivariate random process to be generated. Default is \code{NULL}. 
#' @param ... further arguments for \code{FUN}
#' 
#' @return a matrix or a data frame object
#' 
#' @seealso \code{\link{getVARmodel}}
#' 
#' @title generate
#' @name generate
#' @rdname generate
#' @method generate varest
#' @S3method generate varest
#' @aliases generate generate.varest 
#' 
#' @export
#' 
#' @import RMAWGEN 
#' @examples 
#' 
#' 
#' 
#' 
#' library(RGENERATE)
#' 
#' set.seed(122)
#' NSTEP <- 1000
#' x <- rnorm(NSTEP)
#' y <- x+rnorm(NSTEP)
#' z <- c(rnorm(1),y[-1]+rnorm(NSTEP-1))
#' df <- data.frame(x=x,y=y,z=z)
#' var <- VAR(df,type="none")
#' gg <- generate(var,n=20) 
#'
#' cov <- cov(gg)
#' 
#' ggg <- generate(FUN=rnorm,n=NSTEP,cov=cov)
#' 
#'  
#' library(RMAWGEN)
#' exogen <- as.data.frame(x+5)
#' gpcavar <- getVARmodel(data=df,suffix=NULL,p=3,n_GPCA_iteration=5,
#' n_GPCA_iteration_residuals=5,exogen=exogen)
#' gpcagg <- generate(gpcavar,n=20,exogen=exogen) 
#' 

generate.varest <- function (x,FUN=rnorm,n=100,names=NULL,noise=NULL,exogen=NULL,xprev=NULL,...)  {

	# x is varest objesct
#	out <- require("vars")
#	if (out==FALSE) {
#		
#		print("The function cannot work correctly and returns FALSE")
#		return(out)
#		
#	}
	
	K <-x$K
	p <- x$p
	
	if (x$type!="none") {
	   print("Warning in generate.varest: VAR model type is different from 'none' ")
	   print("Check VAR Model, this method might not work successfully!!")
	}	
	nexogen <- ncol(x$datamat)-K*(p+1)
	## CHECK exogen var!!! 
	if ((is.null(exogen)) & (nexogen!=0)) {
		print("Error in geberate.varest method: exogen variables (predictors) are needed for VAR multirealization") 
		print("Check VAR and exogen variables!!!")
		stop()
	} else if (!is.null(exogen)) if (nexogen!=ncol(exogen)) {
			print("Error in generate.varest method:  corrected exogen variables (predictors) are needed for VAR multirealization") 
			print("Check VAR and exogen variables!!!")
			stop()
	}

	cvar <- coef(x)
	lcvar <- lapply(X=cvar,FUN=function(x){ x[,"Estimate"]})
	mcvar <- matrix(unlist(lcvar),nrow=length(lcvar),byrow=TRUE)
	rownames(mcvar) <- names(lcvar)
	colnames(mcvar) <- names(lcvar[[1]])
	
	
	
	
	if (!is.null(exogen)) n <- nrow(exogen)
	
	if (is.null(xprev)) xprev <- generate(FUN=FUN,n=p,K=K,names=names,cov=summary(x)$covres,...)
	nxprev <- length(xprev)

	if (is.null(noise)) { 
		noise <- generate(FUN=FUN,n=n,K=K,names=names,cov=summary(x)$covres,...)
	##	print(noise)
	
	}
	# MAKE THE GENERATION!!!! 
	
	if (!is.null(exogen)) {
		if (nrow(noise)!=nrow(exogen)) { 
			print("Warning in geberate.varest method:  exogen realizations (predictors) are not equal to noise realizations") 
			nmin <- min(c(nrow(exogen),nrow(noise)))
			exogen <- exogen[1:nmin,]
			noise <- noise[1:nmin,]
			
		}
	}	
	
	out <- noise
# noise and exogen are trasformed into matrix
	noise <- as.matrix(noise)
	if (!is.null(exogen)) exogen <- as.matrix(exogen)
	
	xprev <- as.vector(t(as.matrix(xprev[nrow(xprev):1,])))
	
	
	for (i in 1:nrow(noise)) {
		

		if (is.null(exogen)) { 
			xprevx <- as.vector(xprev)
		} else { 
			xprevx <- as.vector(xprev)

			xprevx <- (c(as.vector(xprev),as.vector((exogen[i,]))))
		}

		xout <- as.vector(mcvar %*% as.matrix(xprevx)) + as.vector(noise[i,])

		if (p>1) { 
			xprev <- c(xout,xprev[1:(K*(p-1))])
		} else { 
			xprev[] <- as.vector(xout)
		}
		
		out[i,] <- xout
	
	}	
	
	return(out)
	
}	

NULL

#' 
#' @rdname generate
#' @method generate varest2
#' @S3method generate varest2
#' @aliases generate 
#' @export

generate.varest2 <- function (x,FUN=rnorm,n=100,names=NULL,noise=NULL,exogen=NULL,xprev=NULL,...) { 

	out <- generate(x=x@VAR,FUN=FUN,n=n,names=names,noise=noise,exogen=exogen,xprev=xprev,...)
	
	return(out)
	
 
}	


NULL

#' 
#' 
#' @param extremes,type  see \code{\link{inv_GPCA}}
#' 
#' @rdname generate
#' @method generate GPCAvarest2
#' @S3method generate GPCAvarest2
#' @aliases generate 
#' @export

generate.GPCAvarest2 <- function (x,FUN=rnorm,n=100,names=NULL,noise=NULL,exogen=NULL,xprev=NULL,extremes=TRUE,type=3,...) { 
	
	# vedere codice RMAWGEN 
	
	K <-x@VAR$K
	p <- x@VAR$p
	
	if (!is.null(exogen)) n <- nrow(exogen)
	if (!is.null(noise))  n <- min(c(nrow(noise),nrow(exogen)))
	
	if (!is.null(noise)) {
		
		cov <- cov(residuals(x))
		noise <- generate(FUN=FUN,n=n,K=K,names=names,cov=cov,...)
		
	}
	
	varresiduals <- inv_GPCA(x=noise,GPCA_param=x@GPCA_residuals,type=type,extremes=extremes)
	out <- generate(x=x@VAR,FUN=FUN,n=n,names=names,noise=varresiduals,exogen=exogen,xprev=xprev,...)
	out <- inv_GPCA(x=out,GPCA_param=x@GPCA_data,type=type,extremes=extremes)
	
	names(out) <- names
	
	return(out)
	
	
}	




#####################################################################
#####################################################################
#####################################################################




# DA VEDERE QUI!!!!!!	
#	
#	K <-var@VAR$K
#	p <- var@VAR$p
#	
#	nexogen <- ncol(var@VAR$datamat)-K*(p+1)
#	if (is.null(noise)) {
#		noise <- array(rnorm(K*nrealization),c(nrealization,K))
#		
#		
#		if (class(var)=="GPCAvarest2") {
#			
#			noise1 <- as.data.frame(t(B %*% t(as.matrix(noise)))) ## DA CORREGGERE QUI!!!
#			
#			noise <- inv_GPCA(x=noise1,GPCA_param=var@GPCA_residuals,type=type,extremes=extremes)
#			B <- diag(1,ncol(B))
#		} 
#	}	else {
#		
#		B <- diag(1,ncol(B))
#	}
#	
#	
#	if ((is.null(exogen)) & (nexogen!=0)) {
#		print("Error exogen variables (predictors) are needed to new VAR multirealization") 
#		print("Check VAR and exogen variables!!!")
#	} else if (!is.null(exogen)) if (nexogen!=ncol(exogen)) {
#			print("Error corrected exogen variables (predictors) are needed to new VAR multirealization") 
#			print("Check VAR and exogen variables!!!")
#		}#if (nexogen!=ncol(exogen)) names
#	
#	out <- array(NA,c(nrealization,K))
#	
#	
#	x <-xprev[1:K]
#	if (p>1) xprev <- xprev[(K+1):(K*p)]
#	
#	
#	for(i in 1:nrow(out)) {
#		
#		
#		if (p>1) xprev <- c(x,xprev[1:(K*(p-1))]) else xprev <- x
#		
#		if (is.null(exogen)) {
#			exogen_data <- NULL
#		} else {
#			exogen_data <- as.vector(t(exogen[i,]))
#		}
#		
#		x <- NewVAReventRealization(var=var@VAR,xprev=xprev,exogen=exogen_data,noise=as.vector(t(noise[i,])),B=B)
#		out[i,] <- x
#		
#		
#	}
#	
#	if (class(var)=="GPCAvarest2") {
#		
#		out <- inv_GPCA(x=out,GPCA_param=var@GPCA_data,type=type,extremes=extremes)
#		
#	} 		
#	
#	return(as.matrix(out))
#	
#	
#	
	
	
	
#	out <- NULL
#	return(out)
	
#}
