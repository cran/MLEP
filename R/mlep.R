mlep <-
function(pedigree, freq){

	pedigree_list = split(pedigree, pedigree[,1])
	lapply(pedigree_list,function(x){
		max_power = sum(x[,6]!=0)+1
		max_term = sum(sapply(0:max_power,function(i)(max_power+1-i)*(max_power+2-i)/2))
		result =.C("mlep",as.integer(nrow(x)),as.integer(unlist(x[,-1])),as.double(freq),as.integer(max_power),
			as.double(rep(0,max_term)),as.integer(rep(0,max_term)),as.integer(0))
		poly = result[[5]][1:result[[7]]]
		attr(poly,"powers") = result[[6]][1:result[[7]]]
		attr(poly,"max_power") = max_power
		poly
	})

}

