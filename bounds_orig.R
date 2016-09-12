# class: a matrix of zeros and ones with number of rows equal to the number of observations and number of columns equal to the number of classess of observations.  Each column of this matrix indicates which particular observations belong to a class.
# delta_prop: a scalar that governs the proportional change between consecutive weights: $\frac{w_i-w_{i-1}}{w_i}\leq\delta_prop\Biggr|\frac{y_i-y_{i-1}}{y_i}\Biggr|$ for $i = 2, ..., n$
# target_shape:
# theta_shape:
# delta_shape:
# grid_length_shape:

bounds = function(y, alpha, beta, class=NULL, delta_prop=NULL, target_shape=NULL, theta_shape=NULL, delta_shape=NULL, grid_length_shape=NULL, solver) {
	
	library("slam")

	n = length(y)
	cvec = c(y, 0)
	
#! Sum of transformed decision variables equal to 1	
	row_ind_1 = rep(1, n)
	col_ind_1 = 1:n
	vals_1 = rep(1, n)
	row_ind_cur = max(row_ind_1)

	if (is.null(class)) {
#! Alpha constraints for the transformed dec. variables
		row_ind_2 = row_ind_cur+sort(rep(1:n, 2))
		col_ind_2 = rep(NA, length(row_ind_2)) 
		col_ind_2[seq(1, 2*n, by=2)] = 1:n
		col_ind_2[seq(2, 2*n, by=2)] = rep(n+1, n)
		vals_2 = rep(c(1, -1/alpha), n)	
		row_ind_cur = max(row_ind_2)

#! Beta constraints for the transformed dec. variables	
		row_ind_3 = row_ind_cur+sort(rep(1:n, 2))
		col_ind_3 = rep(NA, length(row_ind_2)) 
		col_ind_3[seq(1, 2*n, by=2)] = 1:n
		col_ind_3[seq(2, 2*n, by=2)] = rep(n+1, n)
		vals_3 = rep(c(1, -1/beta), n)	
		row_ind_cur = max(row_ind_3)
	
#! Proportional change constraints for the transformed dec. variables			
		if (!is.null(delta_prop)) {
			row_ind_4 = row_ind_cur+sort(rep(1:(n-1), 2))
			col_ind_4 = sort(rep(1:n, 2))[-c(1, 2*n)]
			vals_4 = c(1-delta_prop*(y[2:n]-y[1:(n-1)])/y[2:n], rep(-1, n-1))[rep(c(1, n), n-1)+sort(rep(0:(n-2), 2))]
			row_ind_cur = max(row_ind_4)
			row_ind_5 = row_ind_cur+sort(rep(1:(n-1), 2))
			col_ind_5 = sort(rep(1:n, 2))[-c(1, 2*n)]
			vals_5 = c(1+delta_prop*(y[2:n]-y[1:(n-1)])/y[2:n], rep(-1, n-1))[rep(c(1, n), n-1)+sort(rep(0:(n-2), 2))]
			row_ind_cur = max(row_ind_5)
		}
	}
	
#! Class constraints
	if (!is.null(class)) {
		n_class = ncol(class)
		row_ind_6 = NA
		col_ind_6 = NA
		vals_6 = NA
		for (i in 1:n_class) {
			y_aux = rep(0, n)
			y_aux[which(class[, i]==1)] = y[which(class[, i]==1)] 
			aux = auxiliary(y_aux, alpha[i], beta[i], row_ind_cur) 			
			row_ind_6 = c(row_ind_6, aux$row_ind)
			col_ind_6 = c(col_ind_6, aux$col_ind)
			vals_6 = c(vals_6, aux$vals)
			row_ind_cur = aux$row_ind_cur
		}
		row_ind_6 = row_ind_6[-1]
		col_ind_6 = col_ind_6[-1]
		vals_6 = vals_6[-1]
		row_ind_cur = max(row_ind_6)
	}	
	
#! Shape constraints		
	if (!is.null(delta_shape)) {
		grid = round(seq(n/grid_length_shape, n-n/grid_length_shape, by=n/grid_length_shape))
		col_ind_7 = NA
		row_ind_7 = row_ind_cur+rep(1:length(grid), grid)
		for (i in 1:length(grid)) {
			col_ind_7 = c(col_ind_7, 1:grid[i])
		}
		col_ind_7 = col_ind_7[-1]
		vals_7 = rep(1, sum(grid))
vals_7[cumsum(grid)] = 1/2
		row_ind_cur = max(row_ind_7)
		col_ind_8 = NA
		row_ind_8 = row_ind_cur+rep(1:length(grid), grid)
		for (i in 1:length(grid)) {
			col_ind_8 = c(col_ind_8, 1:grid[i])
		}
		col_ind_8 = col_ind_8[-1]
		vals_8 = rep(1, sum(grid))
vals_8[cumsum(grid)] = 1/2
		row_ind_cur = max(row_ind_8)
	}	
		
#! Build the constraint matrix
	if (is.null(class)) {
		row_ind = c(row_ind_1, row_ind_2, row_ind_3)
		col_ind = c(col_ind_1, col_ind_2, col_ind_3)
		vals = c(vals_1, vals_2, vals_3)
	}
	if (!is.null(delta_prop)) {
		row_ind = c(row_ind_1, row_ind_2, row_ind_3, row_ind_4, row_ind_5)
		col_ind = c(col_ind_1, col_ind_2, col_ind_3, col_ind_4, col_ind_5)
		vals = c(vals_1, vals_2, vals_3, vals_4, vals_5)		
	}	
	if (!is.null(class)) {
		row_ind = c(row_ind_1, row_ind_6)
		col_ind = c(col_ind_1, col_ind_6)
		vals = c(vals_1, vals_6)		
	}
	if (!is.null(delta_shape)) {
		row_ind = c(row_ind_1, row_ind_2, row_ind_3, row_ind_7, row_ind_8)
		col_ind = c(col_ind_1, col_ind_2, col_ind_3, col_ind_7, col_ind_8)
		vals = c(vals_1, vals_2, vals_3, vals_7, vals_8)
	}	
	aux = cbind(row_ind, col_ind, vals)[order(col_ind), ] 
	Amat = simple_triplet_matrix(i = aux[, 1], j = aux[, 2], v = aux[, 3])

#! Define other parameters	
	bvec = c(1, rep(0, 2*n))
	if (!is.null(delta_prop)) {
		bvec = c(1, rep(0, 2*n), rep(0, n-1), rep(0, n-1))		
	}
	if (!is.null(class)) {
		bvec = c(1, rep(0, length(table(row_ind_6))))		
	}
	if (!is.null(delta_shape)) {
		if (target_shape=="normal") {
			bvec_aux_h = NA
			bvec_aux_l = NA
			y_aux = sort(y)
			for (i in grid) {
				bvec_aux_h = c(bvec_aux_h, delta_shape+pnorm(y_aux[i], theta_shape[1], theta_shape[2]))
				bvec_aux_l = c(bvec_aux_l, -delta_shape+pnorm(y_aux[i], theta_shape[1], theta_shape[2]))
			}
			bvec_aux_h = bvec_aux_h[-1]
			bvec_aux_l = bvec_aux_l[-1]
			bvec = c(1, rep(0, 2*n), bvec_aux_h, bvec_aux_l)	
		}
	}
	lb = c(rep(0, n), 0) 
	ub = rep(Inf, n+1)
	sense = c("E", rep("L", n), rep("G", n))
	if (!is.null(delta_prop)) {
		sense = c("E", rep("L", n), rep("G", n), rep("L", n-1), rep("G", n-1))		
	}
	if (!is.null(class)) {
		sense = c("E") 
		for (i in 1:ncol(class)) {
			sense = c(sense, rep("L", length(which(class[, i]!=0))), rep("G", length(which(class[, i]!=0))))
		}
	}
	if (!is.null(delta_shape)) {
		sense = c("E", rep("L", n), rep("G", n), rep("L", length(bvec_aux_h)), rep("G", length(bvec_aux_h)))		
	}
	var_type = rep("C", n+1)	

#! Solve using CPLEX
	if (solver == "cplex") {
		library("Rcplex")
#! Upper bound	
		out = Rcplex(objsense = "max", cvec, Amat, bvec, lb = lb, ub = ub, sense = sense, vtype = var_type, control = list(trace = 0))
		status_h = out$status
		if (status_h == 1) {
			obj_total_h = out$obj
			weights_h = ((out$xopt)[1:n])/((out$xopt)[n+1])
			prices_h = out$extra$lambda
		}
		if (status_h != 1) {
			cat("Optimal solution not found", "\n")
			obj_total_h = NA
			weights_h = NA
			prices_h = NA
		}
#! Lower bound	
		out = Rcplex(objsense = "min", cvec, Amat, bvec, lb = lb, ub = ub, sense = sense, vtype = var_type, control = list(trace = 0))
		status_l = out$status	
		obj_total_l = out$obj
		weights_l = ((out$xopt)[1:n])/((out$xopt)[n+1])
		prices_l = out$extra$lambda
		status_l = out$status
		if (status_l == 1) {
			obj_total_l = out$obj
			weights_l = ((out$xopt)[1:n])/((out$xopt)[n+1])
			prices_l = out$extra$lambda
		}
		if (status_l != 1) {
			cat("Optimal solution not found", "\n")
			obj_total_l = NA
			weights_l = NA
			prices_l = NA
		}
	}
#! Solve using GLPK
	if (solver == "glpk") {
    	library("Rglpk")
    	dir = rep(NA, length(sense))
    	dir[sense=="E"] = '=='
    	dir[sense=="L"] = '<='
    	dir[sense=="G"] = '>='
    	bound = list(lower = list(ind=c(1:length(ub)), val=rep(0,length(ub))), upper = list(ind=c(1:length(ub)), val=ub))
#! Upper bound	
		out = Rglpk_solve_LP(cvec, Amat, dir, bvec, bounds = bound, types = var_type, max = TRUE)
		status_h = out$status
		if (status_h == 0) {
			obj_total_h = out$optimum
			weights_h = ((out$solution)[1:n])/((out$solution)[n+1])
			prices_h = out$extra$lambda
		}
		if (status_h != 0) {
			cat("Optimal solution not found", "\n")		
			obj_total_h = NA
			weights_h = NA
			prices_h = NA
		}		
#! Lower bound	
		out = Rglpk_solve_LP(cvec, Amat, dir, bvec, bounds = bound, types = var_type, max = FALSE)
		status_l = out$status
		if (status_l == 0) {	
			obj_total_l = out$optimum
			weights_l = ((out$solution)[1:n])/((out$solution)[n+1])
			prices_l = out$extra$lambda
		}
		if (status_l != 0) {
			cat("Optimal solution not found", "\n")		
			obj_total_l = NA
			weights_l = NA
			prices_l = NA
		}			
	}	
#! Solve using Gurobi
	if (solver == "gurobi") {
		library("gurobi")
#! Upper bound		
		model = list()
		model$modelsense = 'max'
		model$obj = cvec
		model$A = Amat
		model$sense = rep(NA, length(sense))
		model$sense[sense=="E"] = '='
		model$sense[sense=="L"] = '<='
		model$sense[sense=="G"] = '>='
		model$rhs = bvec
		model$vtype = var_type
		model$ub = ub
		out = gurobi(model)
		status_h = out$status	
		if (status_h == 2) {
			obj_total_h = out$objval
			weights_h = ((out$x)[1:n])/((out$x)[n+1])
			prices_h = out$pi
		}
		if (status_h != 2) {
			cat("Optimal solution not found", "\n")		
			obj_total_h = NA
			weights_h = NA
			prices_h = NA
		}
#! Lower bound			
		model$modelsense = 'min'
		out = gurobi(model)	
		status_l = out$status	
		obj_total_l = NA
		weights_l = NA
		prices_l = NA
		if (status_l == 2) {	
			obj_total_l = out$objval
			weights_l = ((out$x)[1:n])/((out$x)[n+1])
			prices_l = out$pi
		}
		if (status_l != 2) {
			cat("Optimal solution not found", "\n")		
			obj_total_l = NA
			weights_l = NA
			prices_l = NA
		}					
	}			
	
	return(list(status_h = status_h, status_l = status_l, mu_h = obj_total_h, mu_l = obj_total_l, weights_h = weights_h, prices_h = prices_h, weights_l = weights_l, prices_l = prices_l))
}
