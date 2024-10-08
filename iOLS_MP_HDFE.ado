* 15/12/2021 : corrected "cross" in S.E. inversion to increase speed. Note: this required deleting the diagonalization step.
* 15/12/2021 : corrected iteration logical specification
* 16/12/2021 : corrected absorb from varlist to string, as in ppmlhdfe
* 22/12/2021 : coded with matrix multiplication instead of pre-canned program
* 22/12/2021 : added convergence control (limit and maximum)
* 04/01/2022 : added constant + checks for convergence
* 01/02/2022 : drop preserve + post estimation variables
* 04/02/2022 : warm starting point + Correia, Zylkin and Guimarares singleton check
* 20/4/2022 : quiet collinearity + SHOW
* 22/05/2024 : added options to fix delta, rescaled outcome variables, discovered that ln() and log() have different precision levels, changed the iOLS transformation , allowed for noabsorb and absorb in the same package, added convergence checks, added offset, allowed for  increase precision in HDFE as in PPMLHDFE, changed "starting value" of HDFE calls, changed parameter evolution norm 
* 25/05/2024 : create mata functions to increase speed
mata: mata set matacache 5000
mata: mata set matafavor speed
mata: mata set matastrict off
cap program drop iOLS_MP_HDFE
program define iOLS_MP_HDFE, eclass 
syntax varlist [if] [in] [aweight pweight fweight iweight] [, DELta(real 1) LIMit(real 1e-3) OFFset(string) from(name) checkzero(real 1) MAXimum(real 10000) ABSorb(string) SHOW IP FIXED Robust CLuster(string)]        
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


// /*         NORMALIZE DEPENDENT VARIABLE       */
// // Improves convergence in practice
// // Beware: potential problems if y is too small.
// if "`offset'" !="" {
// 	tempvar depvar3
// 	tempvar depvar2
// 	qui : gen `depvar3' = `depvar' / exp(`offset')
// 	quietly : su `depvar3'
// 	local max_y = r(max)
// 	qui : gen `depvar2' = `depvar3' / `max_y'
// }
// else{
// 	tempvar depvar2
// 	quietly : su `depvar'
// 	local max_y = 1
// 	qui : gen `depvar2' = `depvar' / `max_y'
// }
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
quietly: reghdfe `u' `_rhs' [fw=`w']  if `touse' , resid noabsorb  
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
	mata : invXX = invsym(cross(X,X)) 
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
/*         iOLS LOOP       */	
mata: loop_function_nofe(y,X,beta_initial,delta,invXX,criteria,xb_hat,y_tilde,beta_new,past_criteria)
/*         COVARIANCE MATRIX CALCULATION       */	
	mata: alpha = ln(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1])))
	mata : beta_initial[(cols(X)),1] = alpha
	mata: xb_hat = X*beta_initial
	//mata: y_tilde = log(y + delta*exp(xb_hat)) :- (log(delta :+ y:*exp(-xb_hat)) :- ((y:*exp(-xb_hat) :- 1):/(1:+delta)))
	mata: ui = y:*exp(-xb_hat)
	mata: weight = ui:/(1 :+ delta)
	mata: st_numscalar("delta", delta)
	mata: st_local("delta", strofreal(delta))
	cap drop y_tild
	mata: st_store(., st_addvar("double", "y_tild"), "`touse'", y_tilde)
//	quietly: replace y_tild = y_tild + ln(`max_y')
	quietly: reg y_tild `var_list' [`weight'`exp'] if `touse', `option'
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
*********************// ESTIMATION WITH FIXED EFFECTS //************************
********************************************************************************


if "`absorb'" != "" & "`ip'"==""{	
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
quietly: reghdfe `u' `_rhs' [fw=`w']  if `touse' , absorb(`absorb') resid(`e')
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
// di "now"
// cap drop  *_dta*
//  di "`is_cache'"
// local is_cache 
// local cache `'
// local opt "123"
// local opt "gel"
// local usecache 1
// local 1 `is_cache'
// di "`is_cache'"
// local is_cache 1
// local reghdfe_cache 1
// di "`is_cache'"
// local keepids ""
// local savecache `'
// local 1 `usecache'
// di "`usecache'"
// di "`reghdfe_cache'"
// di "`keys'"
// local cachevars ""
// local savecache  0
// di "`savecache'"
// local varlist_cache ""
// local cache ""
// local savecache ""
// local usecache ""
// local cache ""
quietly hdfe `var_list' if `touse'  , absorb(`absorb') generate(M0_)  
local df_a = e(df_a)
cap drop y_tild
quietly gen y_tild = ln(1+`depvar') if `touse'
quietly	hdfe y_tild if `touse'  , absorb(`absorb') generate(Y0_)   tolerance(1e-3) 
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
	mata: xb_hat_N = .
	mata: diff = .
	mata: fe = .
	mata : beta_new = .
	mata past_criteria = .
	local k = 1
	local eps = 1000	
	local almost_conv = 1e-3
/*          LOOP      */	
	mata: loop_function_fe("`touse'", y,xb_hat,xb_hat_M,PX,beta_initial,xb_hat_N,X,diff,Py_tilde,fe,y_tilde,delta,invPXPX,beta_new,criteria,past_criteria)

/*         VARIANCE COVARIANCE CALCULATIONS      */	
	mata: xb_hat_M = PX*beta_initial 
	mata: xb_hat_N = X*beta_initial
	mata: fe = y_tilde - Py_tilde + xb_hat_M - xb_hat_N
	mata: xb_hat = xb_hat_N + fe 
	mata: alpha = ln(mean(y:*exp(-xb_hat)))
	mata: ui = y:*exp(-xb_hat :- alpha)
	mata: weight = ui:/( 1 :+ delta)
  	foreach var in `var_list' {     // rename variables for last ols
	quietly	rename `var' TEMP_`var'
	quietly	rename M0_`var' `var'
	}
cap _crcslbl Y0_ `depvar'
quietly: reg Y0_ `var_list'  if `touse' [`weight'`exp'], `option' noconstant 
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
cap drop iOLS_MP_HDFE_xb_hat
cap drop iOLS_MP_HDFE_fe
cap drop iOLS_MP_HDFE_error
cap drop _reghdfe*
cap drop y_tild
//	mata: xb_hat = (xb_hat :+ alpha :+ ln(`max_y')) // normalization
//	mata: ui = (y*`max_y'):*exp(-xb_hat)
	mata: st_store(., st_addvar("double", "iOLS_MP_HDFE_fe"), "_COPY", fe)
    mata: st_store(., st_addvar("double", "iOLS_MP_HDFE_error"), "_COPY", ui)
   	mata: st_store(., st_addvar("double", "iOLS_MP_HDFE_xb_hat"),"_COPY", xb_hat:+ alpha)
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



********************************************************************************
**************// ESTIMATION WITH INCIDENTAL PARAMETER //************************
********************************************************************************


if "`absorb'" != "" & "`ip'"!=""{
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
quietly: reghdfe `u' `_rhs' [fw=`w']  if `touse' , absorb(`absorb') resid(`e')
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
quietly hdfe `var_list' if `touse' [`weight'] , absorb(`absorb') generate(M0_) 
local df_a = e(df_a)
cap drop y_tild
quietly gen y_tild = ln(1+`depvar') if `touse'
quietly	hdfe y_tild if `touse'  , absorb(`absorb') generate(Y0_)  tolerance(1e-3)  acceleration(sd)  
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
	mata : y_tild = y_tild - mean(y_tild)
}
else {
	mata : beta_initial =invPXPX*cross(PX,Py_tilde)
}
/*         PREPARE LOOP      */	
	mata : criteria = 1000
	mata : xb_hat = .
	mata:  xb_hat_M = .
	mata: xb_hat_N = .
	mata: diff = .
	mata: fe = .
	mata : beta_new = .
	mata past_criteria = .
	local k = 1
	local eps = 1000	
	local almost_conv = 1e-3
	 mata: loop_function_ip("`touse'", y,xb_hat,xb_hat_M,PX,beta_initial,xb_hat_N,X,diff,Py_tilde,fe,y_tilde,delta,invPXPX,beta_new,criteria,past_criteria)
/*         VARIANCE COVARIANCE CALCULATIONS      */	
	mata: xb_hat = X*beta_initial :+ ln(mean(exp(-X*beta_initial):*y))
	mata: ui = y:*exp(-xb_hat)
	mata: weight = ui:/( 1 :+ delta)
  	foreach var in `var_list' {     // rename variables for last ols
	quietly	rename `var' TEMP_`var'
	quietly	rename M0_`var' `var'
	}
cap _crcslbl Y0_ `depvar'
quietly: reg Y0_ `var_list'  if `touse' [`weight'`exp'], `option'  noconstant
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
cap drop iOLS_MP_HDFE_xb_hat
cap drop iOLS_MP_HDFE_fe
cap drop iOLS_MP_HDFE_error
cap drop _reghdfe*
cap drop y_tild
//	mata: xb_hat = (xb_hat :+ alpha :+ ln(`max_y')) // normalization
//	mata: ui = (y*`max_y'):*exp(-xb_hat)
    mata: st_store(., st_addvar("double", "iOLS_MP_HDFE_error"), "_COPY", ui)
   	mata: st_store(., st_addvar("double", "iOLS_MP_HDFE_xb_hat"),"_COPY", xb_hat)
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
cap: mata: mata drop loop_function_nofe()
cap: mata: mata drop loop_function_fe() 
cap: mata: mata drop loop_function_ip() 

mata:
void function loop_function_nofe(y,X,beta_initial,delta,invXX,criteria,xb_hat,y_tilde,beta_new,past_criteria)
{
max = strtoreal(st_local("maximum"))
lim = strtoreal(st_local("limit"))
show = st_local("show")
	for (i=1; i<=max;i++) {
	beta_initial[(cols(X)),1] = ln(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1])))
	xb_hat = X*beta_initial 
	y_tilde = ((y:*exp(-xb_hat) :- 1):/(1:+delta)) + xb_hat
	beta_new = invXX*cross(X,y_tilde)
	past_criteria = criteria
	criteria = max(abs(beta_new:-beta_initial))
	if (past_criteria<criteria) display("Evidence of non-convergence : increasing delta to:") ;;
	if (past_criteria<criteria) delta = delta*1.1 ;;
	if (past_criteria<criteria) delta ;;
	if (past_criteria<criteria) criteria = past_criteria ;;
	if (past_criteria>criteria) beta_initial = beta_new ;;
	if (i == max) display("Maximum number of iterations hit : results are unreliable.") ;; 
 	if (criteria < lim) i=max+1;; // puts an end to the loop 
	if (show != "") criteria;;
	if (mod(i,10)==0) display("Max. Abs. Deviation:") ;;
	if (mod(i,10)==0) 	criteria  ;;
	}
}
end 
  

mata:
void function loop_function_fe(string scalar touse, y,xb_hat,xb_hat_M,PX,beta_initial,xb_hat_N,X,diff,Py_tilde,fe,y_tilde,delta,invPXPX,beta_new,criteria,past_criteria)
{
 max = strtoreal(st_local("maximum"))
 lim = strtoreal(st_local("limit"))
 show = st_local("show")
	for (i=1; i<=max;i++) {
	xb_hat_M = PX*beta_initial 
	xb_hat_N = X*beta_initial
	diff = y_tilde - Py_tilde
	fe = diff + xb_hat_M - xb_hat_N 
	xb_hat =  xb_hat_N + fe   :+ ln(mean(y:*exp(-xb_hat_N - fe)))
	y_tilde = ((y:*exp(-xb_hat) :- 1):/(1:+delta)) + xb_hat 
	stata("cap drop y_tild")
	st_store(., st_addvar("double", "y_tild"), touse, y_tilde-diff)
	stata("cap drop Y0_")
    stata("quietly: hdfe y_tild if \`touse' , absorb(\`absorb') generate(Y0_)  acceleration(sd)   transform(sym)  ")
	st_view(Py_tilde,.,"Y0_",touse)
	beta_new = invPXPX*cross(PX,Py_tilde)
	past_criteria = criteria
	criteria = max(abs(beta_new:-beta_initial))
	if (past_criteria<criteria) display("Evidence of non-convergence : increasing delta to:") ;;
	if (past_criteria<criteria) delta = delta*1.1 ;;
	if (past_criteria<criteria) delta ;;
	if (past_criteria<criteria) criteria = past_criteria ;;
	if (past_criteria>criteria) beta_initial = beta_new ;;
	if (i == max) display("Maximum number of iterations hit : results are unreliable.") ;; 
 	if (criteria < lim) i=max+1;; // puts an end to the loop 
	if (show != "") criteria;;
	if (mod(i,10)==0) display("Max. Abs. Deviation:") ;;
	if (mod(i,10)==0) 	criteria  ;;
	}
}
end
  

mata:
void function loop_function_ip(string scalar touse, y,xb_hat,xb_hat_M,PX,beta_initial,xb_hat_N,X,diff,Py_tilde,fe,y_tilde,delta,invPXPX,beta_new,criteria,past_criteria)
{
max = strtoreal(st_local("maximum"))
lim = strtoreal(st_local("limit"))
show = st_local("show")
	for (i=1; i<=max;i++) {
	diff = y_tilde - Py_tilde 
	xb_hat = X*beta_initial :+ ln(mean(exp(-X*beta_initial):*y))
	y_tilde = y:*exp(-xb_hat):/(1 :+ delta) + xb_hat  
	stata("cap drop y_tild")
	st_store(., st_addvar("double", "y_tild"), touse, y_tilde-diff)
	stata("cap drop Y0_")
    stata("quietly: hdfe y_tild if \`touse' , absorb(\`absorb') generate(Y0_) tolerance(\`almost_conv')  acceleration(sd)   transform(sym)  ")
	st_view(Py_tilde,.,"Y0_",touse)
	beta_new = invPXPX*cross(PX,Py_tilde)
	past_criteria = criteria
	criteria = max(abs(beta_new:-beta_initial))
	if (past_criteria<criteria) display("Evidence of non-convergence : increasing delta to:") ;;
	if (past_criteria<criteria) delta = delta*1.1 ;;
	if (past_criteria<criteria) delta ;;
	if (past_criteria<criteria) criteria = past_criteria ;;
	if (past_criteria>criteria) beta_initial = beta_new ;;
	if (i == max) display("Maximum number of iterations hit : results are unreliable.") ;; 
 	if (criteria < lim) i=max+1;; // puts an end to the loop 
	if (show != "") criteria;;
	if (mod(i,10)==0) display("Max. Abs. Deviation:") ;;
	if (mod(i,10)==0) 	criteria  ;;
	}
}
end  
