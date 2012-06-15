fr <-
function (para, poly) 
{
	Val = sum(sapply(poly, function(x) {
		.C("eval_fr", 
		penetrance = as.double(para), 
		poly = as.double(x), 
		powers = as.integer(attributes(x)$powers), 
		max_power = as.integer(attributes(x)$max_power), 
		length = as.integer(length(x)), 
		result = as.double(0))$result
	}))
	Val
}


