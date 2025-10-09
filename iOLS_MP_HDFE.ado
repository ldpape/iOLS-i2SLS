mata: mata set matacache 5000
mata: mata set matafavor speed
mata: mata set matastrict off
cap program drop iOLS_MP_HDFE
program define iOLS_MP_HDFE, eclass 
syntax varlist [if] [in] [aweight pweight fweight iweight] [, rho(real 1) delta_path(string) LIMit(real 1e-3) WARM OFFset(string) from(name) checkzero(real 1) aweight(varlist) MAXimum(real 10000) POOLsize(real 250) STDopt(string) ABSorb(string) SHOW  FIXED Robust CLuster(string)]        
*------------------------------------------------------------------------------*
*--------------------------     PARSE TEXT     --------------------------------* 
*------------------------------------------------------------------------------*
	marksample touse
	markout `touse'  `cluster', s     
	if  "`robust'" !="" {
		local opt1  = "`robust' "
	}
	if "`cluster'" !="" {
		local opt2 = "vce(cluster `cluster') "
	}
    if "`absorb'" !="" {
		local opt3 = "absorb(`absorb') "
	}
	local option = "`opt1'`opt2' `stdopt'"	
	local list_var `varlist'
	gettoken dvar list_var : list_var // change 
gettoken _rhs list_var : list_var, p("(")
tempvar depvar 
qui: sum `dvar' if `touse' & `dvar'>0
qui: gen `depvar' = `dvar'/r(sd) if `touse' // normalize outcome variable
if r(min)<0{
	di in red "y<0 detected : probable convergence issues"
}
*------------------------------------------------------------------------------*
*---------------------------     DATA CLEANING     ----------------------------* 
*------------------------------------------------------------------------------*
tempvar sep zvar 
qui: gen `sep' =.
qui: gen `zvar' =.
qui: ppmlhdfe `depvar' `_rhs' if  `touse', absorb(`absorb') tagsep(`sep') zvar(`zvar')
if `r(k_omitted)'>0{
di "Collinearity detected - dropping variables :"
}
local drop_list 
local var_list
local i 1
foreach var of varlist `r(fullvarlist)' {
    // Check if current position is *not* in omitted
    if (strpos(" `r(omitted)' ", " `i' ") == 0) {
        local var_list `var_list' `var'
    }
	else{
	    local drop_list `drop_list' `var'
	}
    local ++i
}
di "`drop_list'"
qui: replace `touse' = (`sep'== 0) if `touse' & missing(`sep')==0 // sep is missing when no zeros
//  *** normalize data 
local scales_list ""
   foreach var of varlist `var_list' {	
	qui:  sum `var' if `touse'	
   scalar a_`var' = r(sd)
   qui: replace `var' = `var'/a_`var'	 if `touse' & a_`var'>0
   local invsd = `=r(sd)'
   local scales_list "`scales_list' `invsd'"
   }
mata:    v = strtoreal(tokens("`scales_list'"))
mata:    S = diag(v)
mata:    st_matrix("S", S)
*------------------------------------------------------------------------------*
*------------------------  GPML: CASE WITHOUT FIXED EFFECTS  ------------------* 
*------------------------------------------------------------------------------*
if "`absorb'" == ""{
// *** drop collinear variables 
	tempvar cste 
	qui: gen `cste' = 1
//     quietly: _rmcoll `_rhs'  if `touse', forcedrop 
// 	local var_list `r(varlist)' 
	
*** prepare variables for MATA
	cap drop y_tild
	quietly gen y_tild = ln(1+`depvar') if `touse'
	mata : X=.
	mata : y_tilde =.
	mata : y =.
	mata : st_view(X,.,"`var_list' `cste'","`touse'")
	mata : st_view(y_tilde,.,"y_tild","`touse'")
	mata : st_view(y,.,"`depvar'","`touse'")
	mata w = .
*** weighted regression 
	if "`aweight'"!=""{
	mata : st_view(w,.,"`aweight'","`touse'")
	di in red "Analytical Weights - each obseration is weighted by : sqrt(w)"
}	
	mata : invXX = invsym(cross(X,X))
	if "`aweight'" != ""{
	mata : invXX = invsym(cross(X,w,X))
	}
*** calculate initial values 
capture	 confirm matrix `from'
if _rc==0 {
	mata : beta_initial = st_matrix("`from'")
	mata : beta_initial = beta_initial'
}
else {
	mata : beta_initial = invXX*cross(X,y_tilde) // reg on log (1 + y)/max(y)
}
*** prepare to launch loops
	mata: xb_hat = .
	mata: beta_new = .
	mata: past_criteria = .
	mata: err = .
	mata : criteria = 10000 
	mata : delta = 1
	local k = 1
	local eps = 1000	
	mata: scale_delta = max(y:*exp(-X*beta_initial :-  ln(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1])))))
	mata: stop_crit = 0
*** launch iOLS_delta and/or iOLS_MP 
if "`warm'" != ""{
	if "`delta_path'" == ""{ // if no delta_path, calculate one
    mata: st_local("s1", strofreal(scale_delta*exp(-4)))
    mata: st_local("s2", strofreal(scale_delta*exp(-2)))
    mata: st_local("s3", strofreal(scale_delta))
    mata: st_local("s4", strofreal(scale_delta*exp(1)))
	local delta_path =  "`s1' `s2' `s3' `s4'"
	di "    "
	di "Estimated ẟ - Path: `delta_path'"
	}
capture noisily mata: loop_function_D_nofe(y,X,beta_initial,delta,invXX,criteria,xb_hat,y_tilde,beta_new,past_criteria,w,scale_delta,stop_crit)	
		local rc = _rc
	if `rc' != 0 {
        display as error "/!\ ERROR /!\"
		display as error "Your data has been rescaled. Reload your data before any new analysis."
        error `rc'
    }
}
	mata : delta = `rho' // provided, or equal to 1 
capture noisily mata: loop_function_nofe(y,X,beta_initial,delta,invXX,criteria,xb_hat,y_tilde,beta_new,past_criteria,w,err)
		local rc = _rc
	if `rc' != 0 {
        display as error "/!\ ERROR /!\"
		display as error "Your data has been rescaled. Reload your data before any new analysis."
        error `rc'
    }
*** Covariance Matrix Calculation	
	mata: alpha = ln(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1])))
	mata: beta_initial[(cols(X)),1] = alpha
	mata: xb_hat = X*beta_initial
	mata: ui = y:*exp(-xb_hat)
	mata: weight = ui:/(1 :+ delta)
	mata: st_numscalar("delta", delta)
	mata: st_local("delta", strofreal(delta))
	cap drop y_tild
	mata: st_store(., st_addvar("double", "y_tild"), "`touse'", y_tilde)
	foreach var of varlist `var_list' {	//rescale 
 	quietly: replace `var' = `var'*a_`var'	 if `touse' & a_`var'>0
 	}
		if "`aweight'"==""{
	quietly: reg y_tild `var_list'  if `touse', `option'
		}
		else {
				quietly: reg y_tild `var_list' [aw = `aweight']  if `touse', `option'
		}
	local dof `e(df_r)'
	matrix beta_final = e(b)
	matrix Sigma = e(V)
	mata : Sigma_hat = st_matrix("Sigma")
	mata : Sigma_0 = (cross(X,X))*Sigma_hat*(cross(X,X))
	mata : invXpIWX = invsym(cross(X, weight, X))
	mata : Sigma_tild = invXpIWX*Sigma_0*invXpIWX
 	mata: st_matrix("Sigma_tild", Sigma_tild)
*** Prepare ereturn to return results 
	local names : colnames beta_final
	local nbvar : word count `names'
	mat rownames Sigma_tild = `names' 
    mat colnames Sigma_tild = `names' 
	cap drop _COPY
	quietly: gen _COPY = `touse'
	qui: sum `touse' if `touse'
	scalar Nobs = e(N)
    ereturn post beta_final Sigma_tild , obs(`=e(N)') depname(`dvar') esample(`touse')  dof(`dof') 
    cap drop iOLS_MP_HDFE_xb_hat
	cap drop iOLS_MP_HDFE_error
	cap drop y_tild
*** Restore Variables 
    *mata : xb_hat = (xb_hat :+ ln(mean_y))
    mata : ui  = y:*exp(-xb_hat)
    mata: st_store(., st_addvar("double", "iOLS_MP_HDFE_error"), "_COPY", ui)
    mata: st_store(., st_addvar("double", "iOLS_MP_HDFE_xb_hat"),"_COPY", xb_hat)
	cap drop _COPY
*** export back to stata 	
ereturn scalar delta = `rho'
ereturn  scalar eps =   `eps'
ereturn  scalar niter =  `k'
ereturn scalar df_r = `dof'
ereturn local cmd "iOLS_MP_HDFE"
ereturn local vcetype `option'
ereturn local absvar `absorb'
display in gr _newline(1) _column(1) "Gamma-PML Estimated by iOLS"
di in gr _col(55) "Number of obs = " in ye %8.0f Nobs
ereturn display	
}


*------------------------------------------------------------------------------*
*------------------------ HDFE: CASE WITH FIXED EFFECTS  ----------------------* 
*------------------------------------------------------------------------------*

if "`absorb'" != "" {

*** drop collinear variables 	
// cap drop _COPY_cste
// quietly:	gen _COPY_cste = 1
// quietly: _rmcoll `_rhs' _COPY_cste , forcedrop 
// local var_list `r(varlist)'

*** prepare for looping 
cap drop M0_*
cap drop Y0_*
cap drop xb_hat*
// local is_cache 1
*** weighted regression
	if "`aweight'"!=""{
		di in red "Analytical Weights - each obseration is weighted by : sqrt(w)"
		local aw "aw = `aweight'"
	}
quietly hdfe `var_list' if `touse' [`aw'] , absorb(`absorb') generate(M0_)  acceleration(sd)   transform(sym)
local df_a = e(df_a) // used to correct degrees of freedom 
cap drop y_tild
quietly gen y_tild = ln(1+`depvar') if `touse'
quietly	hdfe y_tild if `touse' [`aw'] , absorb(`absorb') generate(Y0_) acceleration(sd)   transform(sym)  tolerance(1e-3) 
	mata : X=.
	mata : PX=.
	mata : y_tilde =.
	mata : Py_tilde =.
	mata : y =.
	mata : st_view(X,.,"`var_list'","`touse'")
	mata : st_view(PX,.,"M0_*","`touse'")
	mata : st_view(y_tilde,.,"y_tild","`touse'")
	mata : st_view(Py_tilde,.,"Y0_","`touse'")
	mata : st_view(y,.,"`depvar'","`touse'")	
	mata : invPXPX = invsym(cross(PX,PX))
*** initial values 	
capture	 confirm matrix `from'
if _rc==0 {
	mata : beta_initial = st_matrix("`from'")
	mata : beta_initial = beta_initial'
}
else {
	mata : beta_initial = invPXPX*cross(PX,Py_tilde) 
}
*** prepare for loop 
	mata : criteria = 1000
	mata : xb_hat = .
	mata:  xb_hat_M = .
	mata:  fe = .
	mata : beta_new = .
	mata : past_criteria = .
	local k = 1
	local eps = 1000	
	local almost_conv = 1e-2
	mata: k = .
	mata: beta_history = .
	mata: beta_contemporary = .
	mata: stop_crit = 0
	mata: xb_hat_N = .
	mata: diff = .
	mata: err = .
	mata: delta = 1
	mata : scale_delta = max(y:*exp(-PX*beta_initial :- ln(mean(y:*exp(-PX*beta_initial))))) 
*** loop with iOLS_delta and/or iOLS_MP
if "`warm'" != ""{
	if "`delta_path'" == ""{ // if no delta_path, calculate one
    mata: st_local("s1", strofreal(scale_delta*exp(-4)))
    mata: st_local("s2", strofreal(scale_delta*exp(-2)))
    mata: st_local("s3", strofreal(scale_delta))
    mata: st_local("s4", strofreal(scale_delta*exp(1)))
	local delta_path =  "`s1' `s2' `s3' `s4'"
	di "    "
	di "Estimated ẟ - Path: `delta_path'"
	}
capture noisily	mata: loop_function_D("`touse'", y,xb_hat,xb_hat_M,PX,beta_initial,xb_hat_N,X,diff,Py_tilde,fe,y_tilde,delta,invPXPX,beta_new,criteria,past_criteria,scale_delta)
		local rc = _rc
	if `rc' != 0 {
        display as error "/!\ ERROR /!\"
		display as error "Your data has been rescaled. Reload your data before any new analysis."
        error `rc'
    }
}
	mata : delta = `rho'
capture noisily	mata: loop_function_D_fe("`touse'", y,xb_hat,xb_hat_M,PX,beta_initial,xb_hat_N,X,diff,Py_tilde,err,y_tilde,delta,invPXPX,beta_new,criteria,past_criteria)
		local rc = _rc
	if `rc' != 0 {
        display as error "/!\ ERROR /!\"
		display as error "Your data has been rescaled. Reload your data before any new analysis."
        error `rc'
    }
*** variance covariance calculation
 	mata: ui = y:*exp(-xb_hat_M :- log(mean( y:*exp(-xb_hat_M  ))))
	mata: weight =  ui :/ (1 :+ delta)
  	foreach var in `var_list' {     // rename variables for last ols
	quietly	rename `var' TEMP_`var'
	quietly	rename M0_`var' `var'
	}
cap _crcslbl Y0_ `depvar'
	if "`aweight'" == ""{
		quietly: reg Y0_ `var_list'  if `touse' , `option' noconstant 
	}
	else {
		quietly: reg Y0_ `var_list' [aw = `aweight']  if `touse' , `option' noconstant 		
	}
local df_r = e(df_r) - `df_a'
 if "`cluster'" !="" {
 local df_r = e(df_r) 
}
*** report results 	
	matrix b_old = e(b)      
	matrix V_old = e(V)
	matrix beta_final = (inv(S) * b_old')'   // equals diag(s_j) * b_old => rescaling 
	matrix Sigma = (e(df_r) / `df_r')*V_old // correct degrees of freedom 
	foreach var in `var_list' {      // rename variables back
	quietly	rename `var' M0_`var'
	quietly	rename TEMP_`var' `var'
	}
	foreach var of varlist `var_list' {	 // restore scaling 
  	quietly: replace `var' = `var'*a_`var'	 if `touse' & a_`var'>0
  	}
	mata : Sigma_hat = st_matrix("Sigma")
	mata : Sigma_0 = (cross(PX,PX))*Sigma_hat*(cross(PX,PX)) // recover original HAC 
	mata : invXpIWX = invsym(cross(PX, weight,PX)) 
	mata : Sigma_tild = invXpIWX*Sigma_0*invXpIWX
    mata : st_matrix("Sigma_tild", Sigma_tild) // used in practice
	matrix Sigma_tild = inv(S) * Sigma_tild * inv(S)'
	local names : colnames b_old
	local nbvar : word count `names'
	mat rownames Sigma_tild = `names' 
    mat colnames Sigma_tild = `names' 
	mat colnames beta_final = `names' 
	cap drop _COPY
	quietly: gen _COPY = `touse'
	qui: sum `touse' if `touse'
	scalar Nobs = e(N)
    ereturn post beta_final Sigma_tild , obs(`=e(N)') depname(`dvar') esample(`touse')  dof(`df_r')
*** report results (ereturn)
cap drop iOLS_MP_HDFE_error
cap drop _reghdfe*
cap drop y_tild
//mata: ui = y:*exp(-xb_hat_M :- log(mean_y) :- log(mean( y:*exp(-xb_hat_M  )))) // beware, normalization
mata: ui = y:*exp(-xb_hat_M :- log(mean( y:*exp(-xb_hat_M  )))) // beware, normalization
mata: st_store(., st_addvar("double", "iOLS_MP_HDFE_error"), "_COPY", ui)
ereturn scalar rho = `rho'
ereturn scalar eps =   `eps'
ereturn scalar niter =  `k'
ereturn local cmd "iOLS_HDFE_MP"
ereturn local vcetype `option'
ereturn local absvar `absorb'
cap drop _COPY
cap drop Y0_*
cap drop M0_* 
di in gr _col(55) "Number of obs = " in ye %8.0f Nobs
display in gr _newline(1) _column(1) "HDFE iOLS Estimator"
display in gr _newline(1) _column(1) "Absorbed Fixed-Effects: `absorb'"
ereturn display
}
mata: mata drop *
end

*------------------------------------------------------------------------------*
*------------------------ FUNCTIONS REPOSITORY --------------------------------* 
*------------------------------------------------------------------------------*

cap: mata: mata drop loop_function_nofe()
cap: mata: mata drop loop_function_D_nofe()
cap: mata: mata drop loop_function_D() 
cap: mata: mata drop loop_function_D_fe() 

************************** A - NO FIXED EFFECTS ******************************** 
** iOLS_delta loop 
mata:
void function loop_function_D_nofe(y,X,beta_initial,delta,invXX,criteria,xb_hat,y_tilde,beta_new,past_criteria,w,scale_delta,stop_crit)
{
max = strtoreal(st_local("maximum"))
lim = strtoreal(st_local("limit"))
show = st_local("show")
weight = st_local("aweight")
printf("\n")
printf("=========================================================\n")
printf("     Calculating Preliminary Estimate (iOLS-ẟ) \n")
printf("=========================================================\n")
printf("\n")
stop_crit = 0
values = (tokens(st_local("delta_path")))
for (k = 1; k <= length(values); k++) {
    delta = strtoreal(values[k])
 	beta_history = beta_initial
	criteria = 1000
   	printf("Solving for δ equal to: %f\n", delta)
	for (i=1; i<=max;i++) {
	beta_initial[(cols(X)),1] = ln(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1])))
	xb_hat_M = X*beta_initial
	c_hat = mean(log(y :+ delta:*exp(xb_hat_M))) :- mean(xb_hat_M)
	y_tilde = log(y :+ delta:*exp(xb_hat_M)) :- c_hat
	beta_new = invXX*cross(X,y_tilde)
	if (weight!="")  	beta_new = invXX*cross(X, w,  y_tilde) ;;
	if (weight=="")  beta_new = invXX*cross(X,y_tilde) ;;
	criteria = max(abs(beta_new:-beta_initial))
	beta_initial = beta_new
 	if (criteria < 10*lim) i=max+1 ;; // puts an end to the loop 
	if (i == max) display("Maximum number of iterations hit : results are unreliable.") ;; 
	if ((i == 1) & (show != "")) display("Displaying (1) Max. Abs. Dev., (2) Delta, (3) iOLS_delta Step Number") ;; 
	if (show != "") criteria;;
	if (show != "") delta;;
	if (show != "") k;;
	}
beta_contemporary = beta_new 
if (max(abs(beta_contemporary:-beta_history))<lim) k = length(values) + 1
	}
}
end

** iOLS_MP loop 

mata:
void function loop_function_nofe(y,X,beta_initial,delta,invXX,criteria,xb_hat,y_tilde,beta_new,past_criteria,w,err)
{
criteria = 1000
past_criteria = 10000
max = strtoreal(st_local("maximum"))
lim = strtoreal(st_local("limit"))
show = st_local("show")
weight = st_local("aweight")
err = max 
printf("\n")
printf("=========================================================\n")
printf("     Calculating Exact Estimate (iOLS-ρ)\n")
printf("=========================================================\n")
printf("\n")
	for (i=1; i<=max;i++) {
	beta_initial[(cols(X)),1] = ln(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1])))
	xb_hat = X*beta_initial 
	y_tilde = ((y:*exp(-xb_hat) :- 1):/(1:+delta)) + xb_hat
	if (weight=="")  beta_new = invXX*cross(X,y_tilde) ;;
	if (weight!="")  beta_new = invXX*cross(X,w,y_tilde) ;;
	past_criteria = criteria
	criteria = max(abs(beta_new:-beta_initial))
	if (past_criteria<criteria & i>1) delta = delta*2 ;;
	if (past_criteria<criteria & i>1) printf("Convergence Issue - Increasing ρ to: %f\n", delta);; 
	if (past_criteria<criteria & i>1) criteria = past_criteria ;;
	if (past_criteria<criteria & i>1) err = i ;; 
	if (past_criteria>criteria | i==1) beta_initial = beta_new ;;
	if (i == max) display("Maximum number of iterations hit : results are unreliable.") ;; 
 	if (criteria < lim) i=max+1;; // puts an end to the loop 
	if (mod(i,1)==0 & i>1) criteria ;;
	}
if (err<(i-5)) display("Recent Convergence Issue - Results are unreliable") ;; 
printf("\n")
printf("=========================================================\n")
printf("     Final Estimation Results:\n")
printf("=========================================================\n")
printf("\n")
}
end 



************************** B - WITH FIXED EFFECTS ******************************** 
*** iOLS_delta_HDFE loop
mata:
void function loop_function_D(string scalar touse, y,xb_hat,xb_hat_M,PX,beta_initial,xb_hat_N,X,diff,Py_tilde,fe,y_tilde,delta,invPXPX,beta_new,criteria,past_criteria,scale_delta)
{
 max = strtoreal(st_local("maximum"))
 lim = strtoreal(st_local("limit"))
 show = st_local("show")
 weight = st_local("aweight")
printf("\n")
printf("=========================================================\n")
printf("     Calculating Preliminary Estimate (iOLS-ẟ) \n")
printf("=========================================================\n")
printf("\n")
stop_crit = 0 
values = (tokens(st_local("delta_path")))
 for (k = 1; k <= length(values); k++) {
    delta = strtoreal(values[k])
 	beta_history = beta_initial
	criteria = 1000
    printf("Solving for δ equal to: %f\n", delta)
	for (i=1; i<=max ; i++) {
	xb_hat_M = PX*beta_initial
	alpha = log(mean(y:*exp(-xb_hat_M)))
	c_hat = mean(log(y :+ delta:*exp(alpha :+ xb_hat_M))) :- mean(alpha  :+ xb_hat_M)
	diff = y_tilde - Py_tilde 
	y_tilde = log(y :+ delta:*exp(xb_hat_M :+ alpha)) :- c_hat
	stata("cap drop y_tild")
	st_store(., st_addvar("double", "y_tild"), touse, y_tilde-diff)
	stata("cap drop Y0_")
        if (weight=="")  stata("quietly: hdfe y_tild if \`touse' , absorb(\`absorb') generate(Y0_)  acceleration(sd)  tolerance(\`almost_conv')  transform(sym) poolsize(\`poolsize') ") ;;
   	if (weight!="")  stata("quietly: hdfe y_tild if \`touse' [aw = \`aweight'] , absorb(\`absorb') generate(Y0_)  tolerance(\`almost_conv') acceleration(sd)   transform(sym) poolsize(\`poolsize') ") ;;	
	st_view(Py_tilde,.,"Y0_",touse)
	beta_new = invPXPX*cross(PX,Py_tilde)
	criteria = max(abs(beta_new:-beta_initial))
	beta_initial = beta_new
	if (i == max) display("Maximum number of iterations hit : results are unreliable.") ;; 
	if ((i == 1) & (show != "")) display("Displaying (1) Max. Abs. Dev., (2) Delta, (3) iOLS_delta Step Number") ;; 
 	if (criteria < 10*lim) i=max+1 ;; // puts an end to the loop 
	if (show != "") criteria;;
	if (show != "") delta;;
	if (show != "") i;;
								}
beta_contemporary = beta_new 
if (max(abs(beta_contemporary:-beta_history))<lim) k = length(values) + 1
	}
}
end

*** iOLS_MP_HDFE loop

mata:
void function loop_function_D_fe(string scalar touse, y,xb_hat,xb_hat_M,PX,beta_initial,xb_hat_N,X,diff,Py_tilde,err,y_tilde,delta,invPXPX,beta_new,criteria,past_criteria)
{
 criteria = 1000
 past_criteria = 1000
 max = strtoreal(st_local("maximum"))
 lim = strtoreal(st_local("limit"))
 show = st_local("show")
 weight = st_local("aweight")
err = max
printf("\n")
printf("=========================================================\n")
printf("     Calculating Exact Estimate (iOLS-ρ)\n")
printf("=========================================================\n")
printf("\n")
	for (i=1; i<=max;i++) {
 	xb_hat_M = PX*beta_initial 
	alpha = log(mean(y:*exp(-xb_hat_M)))
	diff = y_tilde - Py_tilde
	y_tilde = ((y:*exp(-xb_hat_M :- alpha) :- 1):/(1:+delta)) :+ xb_hat_M  :+ alpha
	stata("cap drop y_tild")
	st_store(., st_addvar("double", "y_tild"), touse, y_tilde-diff)
	stata("cap drop Y0_")
  	if (weight=="")  stata("quietly: hdfe y_tild if \`touse' , absorb(\`absorb') generate(Y0_)  acceleration(sd)   transform(sym) poolsize(\`poolsize')  ") ;;
	if (weight!="")  stata("quietly: hdfe y_tild if \`touse' [aw = \`aweight'] , absorb(\`absorb') generate(Y0_)  acceleration(sd)   transform(sym) poolsize(\`poolsize') ") ;;
	st_view(Py_tilde,.,"Y0_",touse)
	beta_new = invPXPX*cross(PX,Py_tilde)
	past_criteria = criteria
	criteria = max(abs(beta_new:-beta_initial))
 	if (past_criteria<criteria & i>1) delta = delta*2 ;;
	if (past_criteria<criteria & i>1) printf("Convergence Issue - Increasing ρ to: %f\n", delta);; 
 	if (past_criteria<criteria & i>1) criteria = past_criteria ;;
	if (past_criteria<criteria & i>1) err = i ;; 
 	if (past_criteria>criteria | i==1) beta_initial = beta_new ;;
	if (i == max) display("Maximum number of iterations hit : results are unreliable.") ;; 
 	if (criteria < lim) i=max+1;; // puts an end to the loop 
	if (mod(i,1)==0 & i>1) criteria ;;
	}
if (err<(i-5)) display("Recent Convergence Issue - Results are unreliable") ;; 
printf("\n")
printf("=========================================================\n")
printf("     Final Estimation Results:\n")
printf("=========================================================\n")
printf("\n")
}
end


/*


** - Example
* Set the number of individuals (N) and time periods (T)
local N = 5000
local T = 2
set seed 1234
* Create a dataset with all combinations of individuals and time periods
clear 
set obs `N'
gen id = _n
expand `T'
bysort id: gen time = _n


* Generate time-specific effects (gamma)
gen gamma = rnormal(-0.5, 0.5)
gen alpha = rnormal(0.5, 0.5)

* Create a time variable with common shocks across individuals
egen gamma_t = mean(gamma), by(time)
egen alpha_i = mean(alpha), by(id)

* Generate independent variables (X1, X2)
gen X1 =  runiform(0, 1) + alpha_i - gamma_t
gen X2 =  rnormal(0, 1) + X1 - alpha_i + gamma_t
gen Z  =  rnormal(0,1) -gamma_t + alpha_i
gen D  =  rnormal(0, 1) - X1 - X2 - alpha_i + gamma_t - 3*Z
* Generate idiosyncratic errors (epsilon)
gen epsilon = runiform(0, 2)
gen y = exp(2*X1 + 2*X2 + 2*D - gamma_t - ln(abs(0.01+alpha_i))*0.1 )*(epsilon>0.45)
* Create a dependent variable (Y) based on a linear model

gen x1 = X1
 xi: iOLS_MP_HDFE y X1 X2 D x1 x1 x1 i.time ,    warm delta_path(1 10 100)
 xi: iOLS_MP_HDFE y X1 X2 D x1 ,  absorb(time id)   warm delta_path(1 10 100)
replace X1 = X1/10 
  xi: iOLS_MP_HDFE y X1 X2 D x1 ,  absorb(time id)   warm delta_path(1 10 100)

*/
