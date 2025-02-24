mata: mata set matacache 5000
mata: mata set matafavor speed
mata: mata set matastrict off
cap program drop iOLS_MP_HDFE
program define iOLS_MP_HDFE, eclass 
syntax varlist [if] [in] [aweight pweight fweight iweight] [, DELta(real 1) LIMit(real 1e-4) OFFset(string) from(name) checkzero(real 1) aweight(varlist) MAXimum(real 10000) ABSorb(string) SHOW  FIXED Robust CLuster(string)]        
/*         PARSE TEXT       */
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
	local option = "`opt1'`opt2'"	
	local list_var `varlist'
	gettoken depvar list_var : list_var
gettoken _rhs list_var : list_var, p("(")
foreach var of varlist `depvar' `_rhs' {   // drop missing observations
quietly replace `touse' = 0 if missing(`var')
}
tempvar  _mean _ones
qui: gen `_ones' = `depvar' == 0
    if "`absorb'" !="" {
foreach var of varlist `absorb' {   // drop missing observations
cap drop `_mean'
qui: gegen `_mean' = mean(`_ones'), by(`var')
qui: replace `touse' = 0 if `_mean' == 1
}
}
foreach var of varlist `_rhs'{
	qui: gdistinct `var'
	if (r(ndistinct)<3) {
cap drop `_mean'
qui: gegen `_mean' = mean(`_ones'), by(`var')
qui: replace `touse' = 0 if `_mean' == 1	
	}
}
/*         ALGORITHM CHOICE       */	
// first present code for case with no fixed effects
// speed gains from the absence of HDFE calls
********************************************************************************
*********************// ESTIMATION WITHOUT FIXED EFFECTS //*********************
********************************************************************************
if "`absorb'" == ""{
/*         CHECK FOR SEPERATION : CORREIRA CODE      */	
loc tol = 1e-5
tempvar u w xb
quietly: gen `u' =  !`depvar' if `touse'
quietly: su `u'  if `touse', mean
loc K = ceil(r(sum) / `tol' ^ 2)
quietly: gen `w' = cond(`depvar', `K', 1)  if `touse'
quietly: sum `w'
if r(mean)!=0{
while 1 {
cap : reghdfe `u' `_rhs' [fw=`w']  if `touse' , resid noabsorb  
if _rc!=0 {
cap: quietly: reghdfe `u' `_rhs' [fw=`w']  if `touse' , resid noabsorb  
}
quietly: predict double `xb'  if `touse', xbd
quietly:	replace `xb' = 0 if (abs(`xb') < `tol')&(`touse')
quietly:	 cou if (`xb' < 0) & (`touse')
	if !r(N) {
		continue, break
	}
quietly:	replace `u' = max(`xb', 0)  if `touse'
quietly:	drop `xb' `w'
}
quietly: replace `touse'  = (`xb' <= 0) // & (`touse')
}
/*         DROP COLLINEAR VARIABLES      */	
	cap drop _COPY_cste
	gen _COPY_cste = 1
    quietly: _rmcoll `_rhs'  if `touse', forcedrop 
	local var_list `r(varlist)' 
/*         PREPARE iOLS       */	
	cap drop y_tild
	quietly gen y_tild = ln(1+`depvar') if `touse'
	mata : X=.
	mata : y_tilde =.
	mata : y =.
	mata : st_view(X,.,"`var_list' _COPY_cste","`touse'")
	mata : st_view(y_tilde,.,"y_tild","`touse'")
	mata : st_view(y,.,"`depvar'","`touse'")
		mata w = .
	if "`aweight'"!=""{
	mata : st_view(w,.,"`aweight'","`touse'")
	di in red "Analytical Weights - each obseration is weighted by : sqrt(w)"
}	
	mata : invXX = invsym(cross(X,X))
	if "`aweight'" != ""{
	mata : invXX = invsym(cross(X,w,X))

	}
/*         SET INITIAL VALUES       */	
capture	 confirm matrix `from'
if _rc==0 {
	mata : beta_initial = st_matrix("`from'")
	mata : beta_initial = beta_initial'
}
else {
	mata : beta_initial = invXX*cross(X,y_tilde) // reg on log (1 + y)/max(y)
}
/*        INTIALIZE LOOP       */	
	mata: xb_hat = .
	mata: beta_new = .
	mata: past_criteria = .
	mata : criteria = 10000 
	mata : delta = `delta'
	local k = 1
	local eps = 1000	
	mata: scale_delta = max(y:*exp(-X*beta_initial :-  ln(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1])))))
	mata: stop_crit = 0
/*         iOLS LOOP       */	
mata: loop_function_D_nofe(y,X,beta_initial,delta,invXX,criteria,xb_hat,y_tilde,beta_new,past_criteria,w,scale_delta,stop_crit)
mata: loop_function_nofe(y,X,beta_initial,delta,invXX,criteria,xb_hat,y_tilde,beta_new,past_criteria,w)
/*         COVARIANCE MATRIX CALCULATION       */	
	mata: alpha = ln(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1])))
	mata : beta_initial[(cols(X)),1] = alpha
	mata: xb_hat = X*beta_initial
	mata: ui = y:*exp(-xb_hat)
	mata: weight = ui:/(1 :+ delta)
	mata: st_numscalar("delta", delta)
	mata: st_local("delta", strofreal(delta))
	cap drop y_tild
	mata: st_store(., st_addvar("double", "y_tild"), "`touse'", y_tilde)
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
/*         REPORT RESULTS       */	
	local names : colnames beta_final
	local nbvar : word count `names'
	mat rownames Sigma_tild = `names' 
    mat colnames Sigma_tild = `names' 
	cap drop _COPY
	quietly: gen _COPY = `touse'
    ereturn post beta_final Sigma_tild , obs(`=e(N)') depname(`depvar') esample(`touse')  dof(`dof') 
    cap drop iOLS_MP_HDFE_xb_hat
	cap drop iOLS_MP_HDFE_error
	cap drop y_tild
/*         RESTORE PRE-NORMALIZED VALUES       */	
//	mata : xb_hat = (xb_hat :+ ln(`max_y'))
    mata: st_store(., st_addvar("double", "iOLS_MP_HDFE_error"), "_COPY", ui)
    mata: st_store(., st_addvar("double", "iOLS_MP_HDFE_xb_hat"),"_COPY", xb_hat)
	cap drop _COPY
/*         EXPORT CONSTANTS       */	
ereturn scalar delta = `delta'
ereturn  scalar eps =   `eps'
ereturn  scalar niter =  `k'
ereturn scalar df_r = `dof'
ereturn local cmd "iOLS_MP_HDFE"
ereturn local vcetype `option'
di in gr _col(55) "Number of obs = " in ye %8.0f e(N)
ereturn display	
}

********************************************************************************
**************// ESTIMATION WITH FIXED EFFECTS //       ************************
********************************************************************************


if "`absorb'" != "" {
// case with covariates to absorb 
if "`nocheck'" == "1"{
/*         CHECK FOR SEPERATION : CORREIRA CODE      */	
loc tol = 1e-5
tempvar u w xb e
quietly: gen `u' =  !`depvar' if `touse'
quietly: su `u'  if `touse', mean
loc K = ceil(r(sum) / `tol' ^ 2)
quietly: gen `w' = cond(`depvar', `K', 1)  if `touse'
quietly: sum `w'
if r(mean)!=0{
while 1 {
cap: quietly: reghdfe `u' `_rhs' [fw=`w']  if `touse' , absorb(`absorb') resid(`e')
if _rc!=0 {
cap: quietly: reghdfe `u' `_rhs' [fw=`w']  if `touse' , absorb(`absorb') resid(`e')
}
quietly: predict double `xb'  if `touse', xbd
quietly: replace `xb' = 0 if (abs(`xb') < `tol')&(`touse')
quietly:	 cou if (`xb' < 0) & (`touse')
	if !r(N) {
		continue, break
	}
quietly:	replace `u' = max(`xb', 0)  if `touse'
quietly:	drop `xb' `w'
}
quietly: replace `touse'  = (`xb' <= 0) // & (`touse')
}
}
/*         DROP COLLINEAR VARIABLES      */	
cap drop _COPY_cste
quietly:	gen _COPY_cste = 1
quietly: _rmcoll `_rhs' _COPY_cste , forcedrop 
local var_list `r(varlist)'
/*         PREPARE iOLS      */	
cap drop M0_*
cap drop Y0_*
cap drop xb_hat*
local is_cache 1
	if "`aweight'"!=""{
		di in red "Analytical Weights - each obseration is weighted by : sqrt(w)"
		local aw "aw = `aweight'"
	}
quietly hdfe `var_list' if `touse' [`aw'] , absorb(`absorb') generate(M0_)  acceleration(sd)   transform(sym)
local df_a = e(df_a)
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
	mata: delta = `delta'
	mata : invPXPX = invsym(cross(PX,PX))
/*         INITIAL VALUES      */	
capture	 confirm matrix `from'
if _rc==0 {
	mata : beta_initial = st_matrix("`from'")
	mata : beta_initial = beta_initial'
}
else {
	mata : beta_initial = invPXPX*cross(PX,Py_tilde) 
}
/*         PREPARE LOOP      */	
	mata : criteria = 1000
	mata : xb_hat = .
	mata:  xb_hat_M = .
	mata:  fe = .
	mata : beta_new = .
	mata : past_criteria = .
	local k = 1
	local eps = 1000	
	local almost_conv = 1e-1
	mata: k = .
	mata: beta_history = .
	mata: beta_contemporary = .
	mata: stop_crit = 0
	mata: xb_hat_N = .
	mata: diff = .
	mata : scale_delta = max(y:*exp(-PX*beta_initial :- ln(mean(y:*exp(-PX*beta_initial)))))
/*          LOOP      */	
	mata: loop_function_D("`touse'", y,xb_hat,xb_hat_M,PX,beta_initial,xb_hat_N,X,diff,Py_tilde,fe,y_tilde,delta,invPXPX,beta_new,criteria,past_criteria,scale_delta)
	mata: loop_function_D_fe("`touse'", y,xb_hat,xb_hat_M,PX,beta_initial,xb_hat_N,X,diff,Py_tilde,fe,y_tilde,delta,invPXPX,beta_new,criteria,past_criteria)
/*         VARIANCE COVARIANCE CALCULATIONS      */	
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
/*         REPORT RESULTS      */	
	matrix beta_final = e(b) // 
	matrix Sigma = (e(df_r) / `df_r')*e(V)
	foreach var in `var_list' {      // rename variables back
	quietly	rename `var' M0_`var'
	quietly	rename TEMP_`var' `var'
	}
	mata : Sigma_hat = st_matrix("Sigma")
	mata : Sigma_0 = (cross(PX,PX))*Sigma_hat*(cross(PX,PX)) // recover original HAC 
	mata : invXpIWX = invsym(cross(PX, weight,PX)) 
	mata : Sigma_tild = invXpIWX*Sigma_0*invXpIWX
    mata : st_matrix("Sigma_tild", Sigma_tild) // used in practice
	local names : colnames beta_final
	local nbvar : word count `names'
	mat rownames Sigma_tild = `names' 
    mat colnames Sigma_tild = `names' 
	cap drop _COPY
	quietly: gen _COPY = `touse'
    ereturn post beta_final Sigma_tild , obs(`=e(N)') depname(`depvar') esample(`touse')  dof(`df_r')
/*         POST ESTIMATION RESULTS      */
cap drop iOLS_MP_HDFE_error
cap drop _reghdfe*
cap drop y_tild
mata: st_store(., st_addvar("double", "iOLS_MP_HDFE_error"), "_COPY", ui)
ereturn scalar delta = `delta'
ereturn scalar eps =   `eps'
ereturn scalar niter =  `k'
ereturn local cmd "iOLS_HDFE_MP"
ereturn local vcetype `option'
di in gr _col(55) "Number of obs = " in ye %8.0f e(N)
cap drop _COPY
cap drop Y0_*
cap drop M0_* 
ereturn display
}
end

// ADD FUNCTIONS //
cap: mata: mata drop loop_function_nofe()
cap: mata: mata drop loop_function_D_nofe()
cap: mata: mata drop loop_function_D() 
cap: mata: mata drop loop_function_D_fe() 

mata:
void function loop_function_D_nofe(y,X,beta_initial,delta,invXX,criteria,xb_hat,y_tilde,beta_new,past_criteria,w,scale_delta,stop_crit)
{
max = strtoreal(st_local("maximum"))
lim = strtoreal(st_local("limit"))
show = st_local("show")
weight = st_local("aweight")
k = 0
delta = 1
stop_crit = 0
 while (stop_crit == 0) {
 	beta_history = beta_initial
	criteria = 1000
	for (i=1; i<=max;i++) {
	beta_initial[(cols(X)),1] = ln(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1])))
	xb_hat_M = X*beta_initial
	c_hat = mean(log(y :+ delta:*exp(xb_hat_M))) :- mean(xb_hat_M)
	y_tilde = log(y :+ delta:*exp(xb_hat_M)) :- c_hat
	beta_new = invXX*cross(X,y_tilde)
	if (weight!="")  	beta_new = invXX*cross(X, w,  y_tilde) ;;
	if (weight=="")  beta_new = invXX*cross(X,y_tilde) ;;
	past_criteria = criteria
	criteria = max(abs(beta_new:-beta_initial))
	beta_initial = beta_new
 	if (criteria < 1e-2) i=max+1 ;; // puts an end to the loop 
	if (i == max) display("Maximum number of iterations hit : results are unreliable.") ;; 
	if ((i == 1) & (show != "")) display("Displaying (1) Max. Abs. Dev., (2) Delta, (3) iOLS_delta Step Number") ;; 
	if (show != "") criteria;;
	if (show != "") delta;;
	if (show != "") k;;
	}
k = k + 1
beta_contemporary = beta_new 
if (k==1) display("------------- Maximum Absolute Deviations -------------") ;;
(max(abs(beta_contemporary:-beta_history)))
stop_crit = (max(abs(beta_contemporary:-beta_history)))<lim
if (stop_crit==0) delta = exp(k):*scale_delta;;
	}
}
end

mata:
void function loop_function_nofe(y,X,beta_initial,delta,invXX,criteria,xb_hat,y_tilde,beta_new,past_criteria,w)
{
max = strtoreal(st_local("maximum"))
lim = strtoreal(st_local("limit"))
show = st_local("show")
weight = st_local("aweight")
	for (i=1; i<=max;i++) {
	beta_initial[(cols(X)),1] = ln(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1])))
	xb_hat = X*beta_initial 
	y_tilde = ((y:*exp(-xb_hat) :- 1):/(1:+delta)) + xb_hat
	if (weight=="")  beta_new = invXX*cross(X,y_tilde) ;;
	if (weight!="")  	beta_new = invXX*cross(X,w,y_tilde) ;;
	past_criteria = criteria
	criteria = max(abs(beta_new:-beta_initial))
	if (past_criteria<criteria) delta = delta*1.1 ;;
	if (past_criteria<criteria) criteria = past_criteria ;;
	if (past_criteria>criteria) beta_initial = beta_new ;;
 	if (i == 1) display("------------- Final Estimation Step -------------") ;; 	
	if (i == max) display("Maximum number of iterations hit : results are unreliable.") ;; 
 	if (criteria < lim) i=max+1;; // puts an end to the loop 
	if (show != "") criteria;;
	}
}
end 
  

mata:
void function loop_function_D_fe(string scalar touse, y,xb_hat,xb_hat_M,PX,beta_initial,xb_hat_N,X,diff,Py_tilde,fe,y_tilde,delta,invPXPX,beta_new,criteria,past_criteria)
{
 max = strtoreal(st_local("maximum"))
 lim = strtoreal(st_local("limit"))
 show = st_local("show")
 weight = st_local("aweight")
	for (i=1; i<=max;i++) {
 	xb_hat_M = PX*beta_initial
	alpha = log(mean(y:*exp(-xb_hat_M)))
	diff = y_tilde - Py_tilde
	y_tilde = ((y:*exp(-xb_hat_M :- alpha) :- 1):/(1:+delta)) :+ xb_hat_M  :+ alpha
	stata("cap drop y_tild")
	st_store(., st_addvar("double", "y_tild"), touse, y_tilde-diff)
	stata("cap drop Y0_")
    if (weight=="")  stata("quietly: hdfe y_tild if \`touse' , absorb(\`absorb') generate(Y0_)  acceleration(sd)   transform(sym)  ") ;;
	if (weight!="")  stata("quietly: hdfe y_tild if \`touse' [aw = \`aweight'] , absorb(\`absorb') generate(Y0_)  acceleration(sd)   transform(sym)  ") ;;
	st_view(Py_tilde,.,"Y0_",touse)
	beta_new = invPXPX*cross(PX,Py_tilde)
	past_criteria = criteria
	criteria = max(abs(beta_new:-beta_initial))
 	if (past_criteria<criteria) delta = delta*1.1 ;;
 	if (past_criteria<criteria) criteria = past_criteria ;;
 	if (past_criteria>criteria) beta_initial = beta_new ;;
 	if (i == 1) display("------------- Final Estimation Step -------------") ;; 
	if (i == max) display("Maximum number of iterations hit : results are unreliable.") ;; 
 	if (criteria < lim) i=max+1;; // puts an end to the loop 
	if (show != "") criteria;;
	}
}
end


mata:
void function loop_function_D(string scalar touse, y,xb_hat,xb_hat_M,PX,beta_initial,xb_hat_N,X,diff,Py_tilde,fe,y_tilde,delta,invPXPX,beta_new,criteria,past_criteria,scale_delta)
{
 max = strtoreal(st_local("maximum"))
 lim = strtoreal(st_local("limit"))
 show = st_local("show")
 weight = st_local("aweight")
 conv = strtoreal(st_local("almost_conv"))
 k = 0
 delta = 1
 stop_crit = 0
 while (stop_crit == 0) {
 	beta_history = beta_initial
	criteria = 1000
	for (i=1; i<=max ; i++) {
	xb_hat_M = PX*beta_initial
	alpha = log(mean(y:*exp(-xb_hat_M)))
	c_hat = mean(log(y :+ delta:*exp(alpha :+ xb_hat_M))) :- mean(alpha  :+ xb_hat_M)
	diff = y_tilde - Py_tilde 
	y_tilde = log(y :+ delta:*exp(xb_hat_M)) :- c_hat
	stata("cap drop y_tild")
	st_store(., st_addvar("double", "y_tild"), touse, y_tilde-diff)
	stata("cap drop Y0_")
    if (weight=="")  stata("quietly: hdfe y_tild if \`touse' , absorb(\`absorb') generate(Y0_)  acceleration(sd)  tolerance(\`almost_conv')  transform(sym)  ") ;;
    if (weight!="")  stata("quietly: hdfe y_tild if \`touse' [aw = \`aweight'] , absorb(\`absorb') generate(Y0_)  tolerance(\`almost_conv') acceleration(sd)   transform(sym)  ") ;;	
	st_view(Py_tilde,.,"Y0_",touse)
	beta_new = invPXPX*cross(PX,Py_tilde)
	past_criteria = criteria
	criteria = max(abs(beta_new:-beta_initial))
	if (past_criteria<criteria) display("Convergence issue : increasing convergence precision.")
   	if (past_criteria<criteria) st_local("almost_conv", strofreal(max(1e-6,conv*0.1))) ;; // avoid problems due to lax convergence
//   	if (past_criteria<criteria) criteria = past_criteria ;;
//   	if (past_criteria>criteria) beta_initial = beta_new ;;
	beta_initial = beta_new
	if (i == max) display("Maximum number of iterations hit : results are unreliable.") ;; 
	if ((i == 1) & (show != "")) display("Displaying (1) Max. Abs. Dev., (2) Delta, (3) iOLS_delta Step Number") ;; 
 	if (criteria < 1e-2) i=max+1 ;; // puts an end to the loop 
	if (show != "") criteria;;
	if (show != "") delta;;
	if (show != "") k;;
		}
k = k + 1
beta_contemporary = beta_new 
if (k==1) display("------------- Maximum Absolute Deviations -------------") ;;
(max(abs(beta_contemporary:-beta_history)))
stop_crit = (max(abs(beta_contemporary:-beta_history)))<lim
if (stop_crit==0) delta = exp(k):*scale_delta ;;
	}
}
end
