mata: mata set matacache 5000
mata: mata set matafavor speed
mata: mata set matastrict off
cap program drop i2SLS_MP_HDFE
program define i2SLS_MP_HDFE, eclass
syntax varlist [if] [in]  [, rho(real 1) delta_path(string)  ABSorb(varlist) OFFset(string) LIMit(real 1e-3) WARM from(name) nocheck(real 1) MAXimum(real 10000) POOLsize(real 250) ENDog(varlist) INSTR(varlist) STDopt(string) SHOW  FIXED Robust CLuster(string) aweight(varlist)  ]              
*------------------------------------------------------------------------------*
*--------------------------     PARSE TEXT     --------------------------------* 
*------------------------------------------------------------------------------*
	cap drop _reg*
	marksample touse
	markout `touse'  `cluster', s   
	if  "`robust'" !="" {
		local opt1  = "`robust'"
	}
	if "`cluster'" !="" {
		local opt2 = "cluster(`cluster') "
	}
	local option = "`opt1'`opt2' `stdopt'"
	local list_var `varlist'
	gettoken dvar list_var : list_var
	gettoken _rhs list_var : list_var, p("(")
tempvar depvar 
qui: gen `depvar' = `dvar'

*------------------------------------------------------------------------------*
*---------------------------     DATA CLEANING     ----------------------------* 
*------------------------------------------------------------------------------*
// *** drop missing observations 
// foreach var of varlist `depvar' `_rhs' `endog' `instr'{
// quietly  replace `touse' = 0 if missing(`var')	
// }

tempvar sep zvar 
qui: gen `sep' =.
qui: gen `zvar' =.
qui: ppmlhdfe `depvar' `_rhs' `endog' if  `touse', absorb(`absorb') tagsep(`sep') zvar(`zvar')   // `instr'
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

local var_list: list var_list - endog // var_list = exogenous variables 
local var_list: list var_list - instr
qui: replace `touse' = (`sep'== 0) if `touse' & missing(`sep')==0 // sep is missing when no zeros

*qui: sum `depvar' if `touse' & `depvar'>0
*mata: mean_y = `r(mean)' // used to normalize variables during calculations 
*qui:replace `depvar' = `depvar'/r(mean)





	
	
*------------------------------------------------------------------------------*
*------------------------  GPML: CASE WITHOUT FIXED EFFECTS  ------------------* 
*------------------------------------------------------------------------------*
if  "`absorb'" ==""{
	tempvar cste 
	qui: gen `cste' = 1
	local exo_list  `cste' `endog' `var_list'  
	local instr_list `cste' `instr' `var_list' 
	local exogenous `var_list'
*** prepare variables for MATA
	cap drop y_tild
	quietly gen y_tild = ln(1+`depvar') if `touse'
	mata : X=.
	mata : Z=.
	mata : y_tilde =.
	mata : y =.
	mata : st_view(X,.,"`exo_list'","`touse'")
	mata : st_view(Z,.,"`instr_list'","`touse'")
	mata : st_view(y_tilde,.,"y_tild","`touse'")
	mata : st_view(y,.,"`depvar'","`touse'")
	mata : w = . 
*** weighted regression
if "`aweight'"!=""{
	mata : st_view(w,.,"`aweight'","`touse'")
	di in red "Analytical Weights - each obseration is weighted by : sqrt(w)"
}	
*** calculate initial values 
	if "`aweight'"==""{
		mata : invPzX = invsym(cross(X,Z)*invsym(cross(Z,Z))*cross(Z,X))*cross(X,Z)*invsym(cross(Z,Z))
	}
	if "`aweight'"!=""{
	mata : invPzX = invsym(cross(X,w,Z)*invsym(cross(Z,w,Z))*cross(Z,w,X))*cross(X,w,Z)*invsym(cross(Z,w,Z))
	}
capture	 confirm matrix `from'
if _rc==0 {
	mata : beta_initial = st_matrix("`from'")
	mata : beta_initial = beta_initial'
}
else {
	mata : beta_initial = invPzX*cross(Z,y_tilde)
}
*** prepare to launch loops
	mata : delta = 1
	mata : xb_hat = .
	mata : past_criteria = .
	mata : beta_new = .
	mata: k = .
	mata: beta_history = .
	mata: beta_contemporary = .
	mata: stop_crit = 0
	mata : alpha = . 
	mata: c_hat = . 
	mata: err = .
	mata: scale_delta = max(y:*exp(-X*beta_initial :-  ln(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1])))))
	local k = 0
	local eps = 1000	
	mata : criteria = 10000
/*         iOLS LOOP       */	
if "`warm'" != "" {
	if "`delta_path'" == ""{ // if no delta_path, calculate one
    mata: st_local("s1", strofreal(scale_delta*exp(-4)))
    mata: st_local("s2", strofreal(scale_delta*exp(-2)))
    mata: st_local("s3", strofreal(scale_delta))
    mata: st_local("s4", strofreal(scale_delta*exp(1)))
	local delta_path =  "`s1' `s2' `s3' `s4'"
	di "    "
	di "Estimated ẟ - Path: `delta_path'"
	}
	mata:  ivloop_function_D_nofe(y,X,Z,beta_initial,delta,invPzX,criteria,xb_hat,y_tilde,beta_new,past_criteria, stop_crit, beta_history, alpha, c_hat, beta_contemporary,scale_delta,k,w)
}
	mata: delta = `rho'
	mata:  ivloop_function_nofe(y,X,Z,beta_initial,delta,invPzX,criteria,xb_hat,y_tilde,beta_new,past_criteria,w,err)
*** Covariance Matrix Calculation	
	mata: beta_initial[(cols(X)),1] = ln(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1])))
	mata: xb_hat = X*beta_initial
	mata : ui = y:*exp(-xb_hat)
	mata: weight = ui:/(1 :+ delta)
	mata: st_numscalar("delta", delta)
	mata: st_local("delta", strofreal(delta))
	cap drop y_tild
	mata: st_store(., st_addvar("double", "y_tild"), "`touse'", y_tilde)
	if "`aweight'"==""{
		quietly: ivreg2 y_tild `exogenous' (`endog' = `instr') if `touse', `option'
	}
	else {
		quietly: ivreg2 y_tild `exogenous' (`endog' = `instr') [aw=`aweight'] if `touse', `option'
	}
	local dof `e(Fdf2)'
	matrix beta_final = e(b) 
	matrix Sigma = e(V)
	mata : Sigma_hat = st_matrix("Sigma")
	mata : Sigma_0 = (cross(X,Z)*invsym(cross(Z,Z))*cross(Z,X):/rows(X))*Sigma_hat*(cross(X,Z)*invsym(cross(Z,Z))*cross(Z,X):/rows(X)) // recover original HAC 
	mata : invXpPzIWX = invsym(0.5:/rows(X)*cross(X,Z)*invsym(cross(Z,Z))*cross(Z,weight,X)+ 0.5:/rows(X)*cross(X,weight,Z)*invsym(cross(Z,Z))*cross(Z,X))
	mata : Sigma_tild = invXpPzIWX*Sigma_0*invXpPzIWX
	mata : Sigma_tild = (Sigma_tild+Sigma_tild'):/2 
    mata: st_matrix("Sigma_tild", Sigma_tild) // used in practice
*** Prepare ereturn to return results 
	local names : colnames beta_final
	local nbvar : word count `names'
	mat rownames Sigma_tild = `names' 
    mat colnames Sigma_tild = `names' 
	cap drop _COPY
	quietly: gen _COPY = `touse'
	qui: sum `touse' if `touse'
	scalar Nobs = e(N)
    ereturn post beta_final Sigma_tild , obs(`e(N)') depname(`dvar') esample(`touse')  dof(`dof') 
*** Restore Variables 
	cap drop i2SLS_MP_HDFE_xb_hat
	cap drop i2SLS_MP_HDFE_error	
    *mata : xb_hat = (xb_hat :+ ln(mean_y))
    mata : ui  = y:*exp(-xb_hat)
    mata: st_store(., st_addvar("double", "i2SLS_MP_HDFE_error"), "_COPY", ui)
    mata: st_store(., st_addvar("double", "i2SLS_MP_HDFE_xb_hat"),"_COPY", xb_hat)
	cap drop y_tild
	cap drop _COPY
*** export back to stata 	
ereturn scalar delta = `rho'
ereturn  scalar eps =   `eps'
ereturn  scalar niter =  `k'
ereturn scalar widstat = e(widstat) 
ereturn scalar df_r = `dof'
ereturn scalar arf = e(arf)
ereturn local cmd "i2SLS_HDFE_MP"
ereturn local vcetype `option'
ereturn local absvar `absorb'
display in gr _newline(1) _column(1) "Endogenous Gamma-PML Estimated by i2LS"
di in gr _col(55) "Number of obs = " in ye %8.0f Nobs
ereturn display	
}

*------------------------------------------------------------------------------*
*------------------------ HDFE: CASE WITH FIXED EFFECTS  ----------------------* 
*------------------------------------------------------------------------------*
if "`absorb'" != "" {
// *** drop collinear variables 	
// 	local exo_list `endog' var_list 
// 	local instr_list `instr' var_list 
*** prepare for looping 
	cap drop Z0_*
	cap drop E0_*
	cap drop M0_*
	cap drop Y0_*
	cap drop xb_hat*
*** weighted regression
	if "`aweight'"!=""{
		di in red "Analytical Weights - each obseration is weighted by : sqrt(w)"
		local aw "aw = `aweight'"
	}
	if "`var_list'"=="" { // case with no X , only FE 
	quietly hdfe `endog' if `touse' [`aw'] , absorb(`absorb') generate(E0_) acceleration(sd)   transform(sym)
	quietly hdfe `instr' if `touse' [`aw'] , absorb(`absorb') generate(Z0_) acceleration(sd)   transform(sym)
	cap drop y_tild
	quietly gen y_tild = ln(1+`depvar') if `touse'
	quietly	hdfe y_tild  if `touse' [`aw'] , absorb(`absorb') generate(Y0_)  tolerance(1e-3)  acceleration(sd)   transform(sym)
	local df_a = e(df_a)
	local dof_hdfe = e(df_a)
	mata : X=.
	mata : PX=.
	mata : PZ=.
	mata : y_tilde =.
	mata : Py_tilde =.
	mata : y =.
	mata : st_view(X,.,"`endog'","`touse'")
	mata : st_view(PX,.,"E0_*","`touse'")
	mata : st_view(PZ,.,"Z0_*","`touse'")
	mata : st_view(y_tilde,.,"y_tild","`touse'")
	mata : st_view(Py_tilde,.,"Y0_","`touse'")
	mata : st_view(y,.,"`depvar'","`touse'")	
	}
	else { // standard case with both X and FE
	quietly hdfe `var_list'  if `touse'  [`aw'] , absorb(`absorb') generate(M0_) acceleration(sd)   transform(sym)
	quietly hdfe `endog'  if `touse'  [`aw'] , absorb(`absorb') generate(E0_) acceleration(sd)   transform(sym)
	quietly hdfe `instr'  if `touse'  [`aw'] , absorb(`absorb') generate(Z0_) acceleration(sd)   transform(sym)
	cap drop y_tild  
	quietly gen y_tild = ln(1+`depvar') if `touse'
	quietly	hdfe y_tild  if `touse'  [`aw'] , absorb(`absorb') generate(Y0_)  tolerance(1e-3)  acceleration(sd)   transform(sym) 
	local df_a = e(df_a)
	local dof_hdfe = e(df_a)
	mata : X=.
	mata : PX=.
	mata : PZ=.
	mata : y_tilde =.
	mata : Py_tilde =.
	mata : y =.
	mata : st_view(X,.,"`endog' `var_list'","`touse'")
	mata : st_view(PX,.,"E0_* M0_*","`touse'")
	mata : st_view(PZ,.,"Z0_* M0_*","`touse'")
	mata : st_view(y_tilde,.,"y_tild","`touse'")
	mata : st_view(Py_tilde,.,"Y0_","`touse'")
	mata : st_view(y,.,"`depvar'","`touse'")	
	}
*** initial values 	
mata : invPzX = invsym(cross(PX,PZ)*invsym(cross(PZ,PZ))*cross(PZ,PX))*cross(PX,PZ)*invsym(cross(PZ,PZ))
capture	 confirm matrix `from'
if _rc==0 {
	mata : beta_initial = st_matrix("`from'")
	mata : beta_initial = beta_initial'
}
else {
	mata : beta_initial = invPzX*cross(PZ,Py_tilde)
}
*** prepare loop 
	mata: criteria = 1000 // needed to initialize
	local k = 1
	local eps = 1000	
	mata : xb_hat = .
	mata : xb_hat_M = .
	mata : xb_hat_N = .
	mata : diff = .
	mata : fe = .
	mata : beta_new = .
	mata : past_criteria = .
	mata: k = .
	mata: beta_history = .
	mata: beta_contemporary = .
	mata: c_hat = .
	mata: delta = 1 
	mata: err = . 
	mata: stop_crit = 0
	mata : scale_delta = max(y:*exp(-PX*beta_initial :- ln(mean(y:*exp(-PX*beta_initial)))))
	local almost_conv = 1e-2
*** loop with iOLS_delta and/or iOLS_MP
if "`warm'" != "" {
		if "`delta_path'" == ""{ // if no delta_path, calculate one
    mata: st_local("s1", strofreal(scale_delta*exp(-4)))
    mata: st_local("s2", strofreal(scale_delta*exp(-2)))
    mata: st_local("s3", strofreal(scale_delta))
    mata: st_local("s4", strofreal(scale_delta*exp(1)))
	local delta_path =  "`s1' `s2' `s3' `s4'"
	di "    "
	di "Estimated ẟ - Path: `delta_path'"
	}
mata: ivloop_function_D("`touse'", y,xb_hat,xb_hat_M,PX,PZ,beta_initial,xb_hat_N,X,diff,Py_tilde,fe,y_tilde,delta,invPzX,beta_new,criteria,past_criteria,beta_history, c_hat, beta_contemporary, stop_crit,scale_delta)
}
mata: delta = `rho'
mata: ivloop_function_D_fe("`touse'", y,xb_hat,xb_hat_M,PX,PZ,beta_initial,xb_hat_N,X,diff,Py_tilde,err,y_tilde,delta,invPzX,beta_new,criteria,past_criteria)
*** variance covariance calculation
	mata: alpha = log(mean(y:*exp(-xb_hat_M)))
	mata: ui = y:*exp(-xb_hat_M :-alpha)
	mata: weight = ui:/(1 :+ delta)
	foreach var in `var_list' {     // rename variables for last ols
	quietly	rename `var' TEMP_`var'
	quietly	rename M0_`var' `var'
	}	
	foreach var in `instr'  {     // rename variables for last ols
	quietly	rename `var' TEMP_`var'
	quietly	rename Z0_`var' `var'
	}
	foreach var in `endog'  {     // rename variables for last ols
	quietly	rename `var' TEMP_`var'
	quietly	rename E0_`var' `var'
	}
cap _crcslbl Y0_ `depvar' // label Y0 correctly
	if "`aweight'" == ""{
quietly: ivreg2 Y0_ `var_list' (`endog' = `instr')  if `touse' , `option' noconstant   // standard case with X and FE 
	}
	else {
quietly: ivreg2 Y0_ `var_list' (`endog' = `instr') [aw = `aweight'] if `touse' , `option' noconstant   // standard case with X and FE 		
	}
local df_r = e(Fdf2) - `df_a'
 if "`cluster'" !="" {
 local df_r = e(Fdf2)
}
	foreach var in `var_list' {      // rename variables back
	quietly	rename `var' M0_`var'
	quietly	rename TEMP_`var' `var'
	}
	foreach var in  `instr' {      // rename variables back
	quietly	rename `var' Z0_`var'
	quietly	rename TEMP_`var' `var'
	}
	foreach var in `endog'{      // rename variables back
	quietly	rename `var' E0_`var'
	quietly	rename TEMP_`var' `var'
	}
	matrix beta_final = e(b) 
	matrix Sigma = (e(Fdf2) / `df_r')*e(V)
	mata : Sigma_hat = st_matrix("Sigma")
	mata : Sigma_0 = (cross(PX,PZ)*invsym(cross(PZ,PZ))*cross(PZ,PX):/rows(PX))*Sigma_hat*(cross(PX,PZ)*invsym(cross(PZ,PZ))*cross(PZ,PX):/rows(PX)) // recover original HAC 
	mata : invXpPzIWX = invsym(0.5:/rows(PX)*cross(PX,PZ)*invsym(cross(PZ,PZ))*cross(PZ,weight,PX)+ 0.5:/rows(PX)*cross(PX,weight,PZ)*invsym(cross(PZ,PZ))*cross(PZ,PX))
	mata : Sigma_tild = invXpPzIWX*Sigma_0*invXpPzIWX
	mata : Sigma_tild = (Sigma_tild+Sigma_tild'):/2 
    mata: st_matrix("Sigma_tild", Sigma_tild) // used in practice
*** report results (ereturn)
	local names : colnames beta_final
	local nbvar : word count `names'
	mat rownames Sigma_tild = `names' 
    mat colnames Sigma_tild = `names' 
	local dof_final = e(df r)- `dof_hdfe'
	cap drop _COPY
	quietly: gen _COPY = `touse'
	qui: sum `touse' if `touse'
	scalar Nobs = e(N)
    ereturn post beta_final Sigma_tild , obs(`=e(N)') depname(`dvar') esample(`touse')  dof(`df_r')
	cap drop i2SLS_MP_HDFE_error
	cap drop _reghdfe*
	cap drop y_tild
*** report normalized variables 
*mata: ui = y:*exp(-xb_hat_M :- log(mean_y) :- log(mean( y:*exp(-xb_hat_M  ))))
mata: ui = y:*exp(-xb_hat_M  :- log(mean( y:*exp(-xb_hat_M  ))))
mata: st_store(., st_addvar("double", "i2SLS_MP_HDFE_error"), "_COPY", ui)
		cap drop _COPY
*** ereturn scalars 
ereturn scalar delta = `rho'
ereturn  scalar eps =   `eps'
ereturn  scalar niter =  `k'
ereturn local cmd "i2SLS_HDFE"
ereturn local vcetype `option'
ereturn local absvar `absorb'
di in gr _col(55) "Number of obs = " in ye %8.0f Nobs
display in gr _newline(1) _column(1) "HDFE i2SLS Estimator"
display in gr _newline(1) _column(4) "Absorbed Fixed-Effects: `absorb'"
ereturn display
**** drop tempvars 
	cap drop E0_*
	cap drop Z0_*
	cap drop M0_* 
	cap drop Y0_*
	cap drop xb_hat*
	}
mata: mata drop *
end


*------------------------------------------------------------------------------*
*------------------------ FUNCTIONS REPOSITORY --------------------------------* 
*------------------------------------------------------------------------------*

cap: mata: mata drop ivloop_function_D_nofe() 
cap: mata: mata drop ivloop_function_nofe()
cap: mata: mata drop ivloop_function_D() 
cap: mata: mata drop ivloop_function_D_fe() 

************************** A - NO FIXED EFFECTS ******************************** 

** i2SLS_delta loop 

mata:
void function ivloop_function_D_nofe(y,X,Z,beta_initial,delta,invPzX,criteria,xb_hat,y_tilde,beta_new,past_criteria, stop_crit, beta_history, alpha, c_hat, beta_contemporary,scale_delta,k,w)
{
max = strtoreal(st_local("maximum"))
lim = strtoreal(st_local("limit"))
show = (st_local("show"))
weight = st_local("aweight")
printf("\n")
printf("=========================================================\n")
printf("     Calculating Preliminary Estimate (i2SLS-ẟ) \n")
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
	beta_initial[(cols(X)),1] = ln(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1])))
	xb_hat = X*beta_initial 
	c_hat = mean(log(y :+ delta:*exp(xb_hat))) :- mean(xb_hat)
	y_tilde = log(y :+ delta:*exp(xb_hat)) :- c_hat
	if (weight=="")  beta_new = invPzX*cross(Z,y_tilde) ;;
	if (weight!="")  beta_new = invPzX*cross(Z, w , y_tilde)  ;; 
	past_criteria = criteria
	criteria = max(abs(beta_new:-beta_initial))
	beta_initial = beta_new
 	if (criteria < 10*lim) i=max+1 ;; // puts an end to the loop 
	if (i == max) display("Maximum number of iterations hit : results are unreliable.") ;; 
	if ((i == 1) & (show != "")) display("Displaying (1) Max. Abs. Dev., (2) Delta, (3) i2SLS_delta Step Number") ;; 
	if (show != "") criteria;;
	if (show != "") delta;;
	if (show != "") k;;
	}
beta_contemporary = beta_new 
if (max(abs(beta_contemporary:-beta_history))<lim) k = length(values) + 1
	}
}
end

** i2SLS_MP loop 

mata:
void function ivloop_function_nofe(y,X,Z,beta_initial,delta,invPzX,criteria,xb_hat,y_tilde,beta_new,past_criteria,w,err)
{
criteria = 1000
past_criteria = 10000
max = strtoreal(st_local("maximum"))
lim = strtoreal(st_local("limit"))
show = (st_local("show"))
weight = st_local("aweight")
err = max 
printf("\n")
printf("=========================================================\n")
printf("     Calculating Exact Estimate (i2SLS-ρ)\n")
printf("=========================================================\n")
printf("\n")
	for (i=1; i<=max;i++) {
	beta_initial[(cols(X)),1] = ln(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1])))
	xb_hat = X*beta_initial
	y_tilde = ((y:*exp(-xb_hat) :- 1):/(1:+delta)) + xb_hat  
	if (weight=="")  beta_new = invPzX*cross(Z,y_tilde) ;;
 	if (weight!="")  beta_new = invPzX*cross(Z, w , y_tilde)  ;; 
	past_criteria = criteria
	criteria = max(abs(beta_new:-beta_initial))
	if (past_criteria<criteria & i>1) delta = delta*2 ;;
	if (past_criteria<criteria & i>1) printf("Convergence Issue - Increasing ρ to: %f\n", delta);; 
	if (past_criteria<criteria & i>1) err = i ;; 
	if (past_criteria<criteria & i>1) criteria = past_criteria ;;
	if (past_criteria>criteria | i==1 ) beta_initial = beta_new ;;
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

*** i2SLS_delta_HDFE loop

mata:
void function ivloop_function_D(string scalar touse, y,xb_hat,xb_hat_M,PX,PZ,beta_initial,xb_hat_N,X,diff,Py_tilde,fe,y_tilde,delta,invPzX,beta_new,criteria,past_criteria,beta_history, c_hat, beta_contemporary, stop_crit,scale_delta)
{
max = strtoreal(st_local("maximum"))
lim = strtoreal(st_local("limit"))
show = (st_local("show"))
weight = st_local("aweight")
conv = strtoreal(st_local("almost_conv"))
printf("\n")
printf("=========================================================\n")
printf("     Calculating Preliminary Estimate (i2SLS-ẟ) \n")
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
	diff = y_tilde - Py_tilde
	alpha = log(mean(y:*exp(-xb_hat_M)))
	c_hat = mean(log(y :+ delta:*exp(alpha :+ xb_hat_M))) :- mean(alpha  :+ xb_hat_M)
	y_tilde = log(y :+ delta:*exp(xb_hat_M :+ alpha)) :- c_hat
	stata("cap drop y_tild")
	st_store(., st_addvar("double", "y_tild"), touse, y_tilde - diff)
	stata("cap drop Y0_")
	if (weight=="") stata("quietly: hdfe y_tild if \`touse' , absorb(\`absorb') generate(Y0_)  tolerance(\`almost_conv')  acceleration(sd)   transform(sym) poolsize(\`poolsize') ")  ;;
        if (weight!="") stata("quietly: hdfe y_tild if \`touse' [aw = \`aweight'] , absorb(\`absorb') generate(Y0_)  tolerance(\`almost_conv')  acceleration(sd)   transform(sym) poolsize(\`poolsize')")  ;;
	st_view(Py_tilde,.,"Y0_",touse)
	beta_new = invPzX*cross(PZ,Py_tilde)
	criteria = max(abs(beta_new:-beta_initial))
	beta_initial = beta_new
	if (i == max) display("Maximum number of iterations hit : results are unreliable.") ;; 
	if ((i == max) & (show != "")) display("Displaying (1) Max. Abs. Dev., (2) Delta, (3) i2SLS_delta Step Number") ;; 
 	if (criteria < 10*lim) i=max+1 ;; // puts an end to the loop 
	if (show != "") criteria;;
	if (show != "") delta;;
	if (show != "") k;;
							}
beta_contemporary = beta_new 
if (max(abs(beta_contemporary:-beta_history))<lim) k = length(values) + 1 ;;
	}
}
end

*** i2SLS_MP_HDFE loop


mata:
void function ivloop_function_D_fe(string scalar touse, y,xb_hat,xb_hat_M,PX,PZ,beta_initial,xb_hat_N,X,diff,Py_tilde,err,y_tilde,delta,invPzX,beta_new,criteria,past_criteria)
{
criteria = 1000
past_criteria = 1000
max = strtoreal(st_local("maximum"))
lim = strtoreal(st_local("limit"))
show = (st_local("show"))
weight = st_local("aweight")
err = max
printf("\n")
printf("=========================================================\n")
printf("     Calculating Exact Estimate (i2SLS-ρ)\n")
printf("=========================================================\n")
printf("\n")
	for (i=1; i<=max ; i++) {
	xb_hat_M = PX*beta_initial 
	alpha = log(mean(y:*exp(-xb_hat_M)))
	diff = y_tilde - Py_tilde
   	y_tilde = ((y:*exp(-xb_hat_M :- alpha)  :- 1):/(1:+delta)) + xb_hat_M :+ alpha 
	stata("cap drop y_tild")
	st_store(., st_addvar("double", "y_tild"), touse, y_tilde-diff)
	stata("cap drop Y0_")
   	 if (weight=="") stata("quietly: hdfe y_tild if \`touse' , absorb(\`absorb') generate(Y0_)    acceleration(sd)   transform(sym) poolsize(\`poolsize') ")  ;;
   	 if (weight!="") stata("quietly: hdfe y_tild if \`touse' [aw = \`aweight'] , absorb(\`absorb') generate(Y0_)   acceleration(sd)   transform(sym) poolsize(\`poolsize')")  ;;
	st_view(Py_tilde,.,"Y0_",touse)
	beta_new = invPzX*cross(PZ,Py_tilde)
	past_criteria = criteria
	criteria = max(abs(beta_new:-beta_initial))
	if (past_criteria<criteria & i>1) delta = delta*2 ;;
	if (past_criteria<criteria & i>1) printf("Convergence Issue - Increasing ρ to: %f\n", delta);; 
	if (past_criteria<criteria & i>1) criteria = past_criteria ;;
	if (past_criteria<criteria & i>1) err = i ;; 
	if (past_criteria>criteria | i==1) beta_initial = beta_new ;;
	if (i == max) display("Maximum number of iterations hit : results are unreliable.") ;; 
 	if (criteria < lim) i=max+1 ;; // puts an end to the loop 
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
local N = 10000
local T = 2
set seed 1234
* Create a dataset with all combinations of individuals and time periods
clear 
set obs `N'
gen id = _n
expand `T'
bysort id: gen time = _n

* Generate individual-specific effects (alpha)
gen alpha = rnormal(0.2, 0.5)

* Generate time-specific effects (gamma)
gen gamma = rnormal(0.2, 0.5)

* Create a time variable with common shocks across individuals
egen gamma_t = mean(gamma), by(time)
egen alpha_i = mean(alpha), by(id)

* Generate independent variables (X1, X2)
gen X1 =  runiform(0, 1) + alpha_i - gamma_t
gen X2 =  rnormal(0, 1) + X1 - alpha_i + gamma_t
gen Z  =  rnormal(0,1) -0.1*gamma_t + 0.1*alpha_i
gen D  =  rnormal(0, 1) + X1 - X2 - alpha_i - gamma_t - Z
* Generate idiosyncratic errors (epsilon)
gen epsilon = runiform(0, 2)
gen Y = ( 0.5*abs(gamma_t) + 0.1*abs(alpha_i))*exp(-X1 + X2 + D)*epsilon
gen wvar = (uniform())*1000 // random weights 
* Create a dependent variable (Y) based on a linear model
 xi: i2SLS_MP_HDFE Y X1   i.time  ,  endog(D) instr(Z) warm 
i2SLS_MP_HDFE Y X1 X2 , absorb(time id)  endog(D) instr(Z)  warm
*/
