grr <-
function(para, poly, minimize=F){

	Val = apply(sapply(poly, function(x){
		.C("eval_grr",            
			penetrance=as.double(para),
			poly = as.double(x),
			powers = as.integer(attributes(x)$powers),
			max_power = as.integer(attributes(x)$max_power),
			length = as.integer(length(x)), 
			result=as.double(rep(0,3)) 
			)$result
		}), 1, sum)
	if(minimize){
		-Val
	}else{
		Val
	}
}

