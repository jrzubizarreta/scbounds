auxiliary = function(y, alpha, beta, row_ind_cur) {
	n = length(which(y!=0))
	
#! Alpha bounds for the transformed dec. variables
	row_ind_a = row_ind_cur+sort(rep(1:n, 2))
	col_ind_a = rep(NA, length(row_ind_a)) 
	col_ind_a[seq(1, 2*n, by=2)] = which(y!=0)
	col_ind_a[seq(2, 2*n, by=2)] = rep(length(y)+1, n)
	vals_a = rep(c(1, -1/alpha), n)	
	row_ind_cur = max(row_ind_a)

#! Beta bounds for the transformed dec. variables	
	row_ind_b = row_ind_cur+sort(rep(1:n, 2))
	col_ind_b = rep(NA, length(row_ind_a)) 
	col_ind_b[seq(1, 2*n, by=2)] = which(y!=0)
	col_ind_b[seq(2, 2*n, by=2)] = rep(length(y)+1, n)
	vals_b = rep(c(1, -1/beta), n)	
	row_ind_cur = max(row_ind_b)
	
	row_ind = c(row_ind_a, row_ind_b)
	col_ind = c(col_ind_a, col_ind_b)
	vals = c(vals_a, vals_b)
	
	return(list(row_ind = row_ind, col_ind = col_ind, vals = vals, row_ind_cur = row_ind_cur))	
}