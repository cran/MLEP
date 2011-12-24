persp_penetrance <-
function(poly, fixed, value, ...){

	if(fixed == "alpha"){
		zval = t(sapply(0:10/10,function(x){
			sapply(0:10/10,function(y){
				ifelse(x<y,NA,fr(c(value,x,y),poly))
			})
		}))
		xlab="beta"
		ylab="gamma"
	}
	else if(fixed == "beta"){
		zval = t(sapply(0:10/10,function(x){
			sapply(0:10/10,function(y){
				ifelse(x<y,NA,fr(c(x,value,y),poly))
			})
		}))
		xlab="alpha"
		ylab="gamma"
	}
	else{
		zval = t(sapply(0:10/10,function(x){
			sapply(0:10/10,function(y){
				ifelse(x<y,NA,fr(c(x,y,value),poly))
			})
		}))
		xlab="alpha"
		ylab="beta"
	}

	zval[zval=="-Inf"]=NA
	persp(0:10/10,0:10/10,zval,
		xlab=xlab,ylab=ylab,zlab="Log Likelihood",ticktype = "detailed",...)

}

