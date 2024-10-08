* 15/12/2021 : corrected "cross" in S.E. which is complicated by the symmetrization
* 15/12/2021 : corrected iteration logical specification
* 16/12/2021 : corrected absorb from varlist to string, as in ppmlhdfe
* 22/12/2021 : coded with matrix multiplication instead of pre-canned program
* 22/12/2021 : added convergence control (limit and maximum)
* 04/01/2022 : added constant + checks for convergence + corrected problem with collinear variables affecting final 2SLS
* 21/01/2022 : added symmetrization + check for singleton / existence using PPML + correction S.E. + syntax change
* 22/01/2022 : apparently, new syntax does not drop missing obs.
*3/2/2022 : drop preserve + add singleton selection based on Correia, Zylkin and Guimaraes.
* 20/4/2022 : quietly collinearity + SHOW option 
* 22/05/2024 : added options to fix delta, rescaled outcome variables, discovered that ln() and log() have different precision levels, changed the iOLS transformation , allowed for noabsorb and absorb in the same package, added convergence checks, added offset, allowed for  increase precision in HDFE as in PPMLHDFE, changed "starting value" of HDFE calls, changed parameter evolution norm 
* 25/05/2024 : create mata functions to increase speed
mata: mata set matacache 5000
mata: mata set matafavor speed
mata: mata set matastrict off
cap program drop i2SLS_MP_HDFE
program define i2SLS_MP_HDFE, eclass
syntax varlist [if] [in] [aweight pweight fweight iweight] [, DELta(real 1)  ABSorb(varlist) OFFset(string) LIMit(real 1e-3) from(name)   nocheck(real 1) MAXimum(real 10000) ENDog(varlist) INSTR(varlist) SHOW IP FIXED Robust CLuster(string)]              
/*         PARSE TEXT       */
	cap drop _reg*
	marksample touse
	markout `touse'  `cluster', s   
	if  "`robust'" !="" {
		local opt1  = "`robust' "
	}
	if "`cluster'" !="" {
		local opt2 = "cluster(`cluster') "
	}
	local option = "`opt1'`opt2'"
	local list_var `varlist'
	gettoken depvar list_var : list_var
	gettoken _rhs list_var : list_var, p("(")
foreach var of varlist `depvar' `_rhs' `endog' `instr'{
quietly  replace `touse' = 0 if missing(`var')	
}
/*         NORMALIZE DEPENDENT VARIABLE       */
// Improves convergence in practice
// Beware: potential problems if y is too small.
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
// 	local max_y = r(max)
// 	qui : gen `depvar2' = `depvar' / `max_y'
// }
/*         ALGORITHM CHOICE       */	
// first present code for case with no fixed effects
// speed gains from the absence of HDFE calls


********************************************************************************
*********************// ESTIMATION WITHOUT FIXED EFFECTS //*********************
********************************************************************************

if  "`absorb'" ==""{
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
	tempvar cste
	gen `cste' = 1
	quietly:   _rmcoll `_rhs' `cste' if `touse', forcedrop 
	local var_list `endog' `r(varlist)' `cste'  
	local instr_list `instr' `r(varlist)' `cste' 
	local exogenous `r(varlist)'
/*         PREPARE iOLS       */
	cap drop y_tild
	quietly gen y_tild = ln(1+`depvar') if `touse'
	mata : X=.
	mata : Z=.
	mata : y_tilde =.
	mata : y =.
	mata : st_view(X,.,"`var_list'","`touse'")
	mata : st_view(Z,.,"`instr_list'","`touse'")
	mata : st_view(y_tilde,.,"y_tild","`touse'")
	mata : st_view(y,.,"`depvar'","`touse'")
	mata : invPzX = invsym(cross(X,Z)*invsym(cross(Z,Z))*cross(Z,X))*cross(X,Z)*invsym(cross(Z,Z))
/*         SET INITIAL VALUES       */	
capture	 confirm matrix `from'
if _rc==0 {
	mata : beta_initial = st_matrix("`from'")
	mata : beta_initial = beta_initial'
}
else {
	mata : beta_initial = invPzX*cross(Z,y_tilde)
}
/*        INTIALIZE LOOP       */	
	mata : delta = `delta'
	mata : xb_hat = .
	mata : past_criteria = .
	mata : beta_new = .
	local k = 1
	local eps = 1000	
	mata : criteria = 10000
/*         iOLS LOOP       */	

	mata:  ivloop_function_nofe(y,X,Z,beta_initial,delta,invPzX,criteria,xb_hat,y_tilde,beta_new,past_criteria)

/*         VARIANCE COVARIANCE CALCULATIONS      */	
	mata: beta_initial[(cols(X)),1] = ln(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1])))
	mata: xb_hat = X*beta_initial
	mata : ui = y:*exp(-xb_hat)
	mata: weight = ui:/(1 :+ delta)
	mata: st_numscalar("delta", delta)
	mata: st_local("delta", strofreal(delta))
	cap drop y_tild
	mata: st_store(., st_addvar("double", "y_tild"), "`touse'", y_tilde)
//	quietly: replace y_tild = y_tild + ln(`max_y')
    quietly: ivreg2 y_tild `exogenous' (`endog' = `instr') [`weight'`exp'] if `touse', `option'
	local dof `e(Fdf2)'
	matrix beta_final = e(b) // 	mata: st_matrix("beta_final", beta_new)
	matrix Sigma = e(V)
	mata : Sigma_hat = st_matrix("Sigma")
	mata : Sigma_0 = (cross(X,Z)*invsym(cross(Z,Z))*cross(Z,X):/rows(X))*Sigma_hat*(cross(X,Z)*invsym(cross(Z,Z))*cross(Z,X):/rows(X)) // recover original HAC 
	mata : invXpPzIWX = invsym(0.5:/rows(X)*cross(X,Z)*invsym(cross(Z,Z))*cross(Z,weight,X)+ 0.5:/rows(X)*cross(X,weight,Z)*invsym(cross(Z,Z))*cross(Z,X))
	mata : Sigma_tild = invXpPzIWX*Sigma_0*invXpPzIWX
	mata : Sigma_tild = (Sigma_tild+Sigma_tild'):/2 
    mata: st_matrix("Sigma_tild", Sigma_tild) // used in practice
/*         REPORT RESULTS       */	
	local names : colnames beta_final
	local nbvar : word count `names'
	mat rownames Sigma_tild = `names' 
    mat colnames Sigma_tild = `names' 
	cap drop _COPY
	quietly: gen _COPY = `touse'
    ereturn post beta_final Sigma_tild , obs(`e(N)') depname(`depvar') esample(`touse')  dof(`dof') 
/*         RESTORE PRE-NORMALIZED VALUES       */	
	//mata : xb_hat = (xb_hat+log(`y_max'))
	//mata : ui = -(y*`y_max'):*exp(-xb_hat)
    mata: st_store(., st_addvar("double", "i2SLS_MP_HDFE_error"), "_COPY", ui)
    mata: st_store(., st_addvar("double", "i2SLS_MP_HDFE_xb_hat"),"_COPY", xb_hat)
	cap drop i2SLS_MP_HDFE_xb_hat
	cap drop i2SLS_MP_HDFE_error
	cap drop y_tild
	cap drop _COPY
/*         EXPORT CONSTANTS       */	
ereturn scalar delta = `delta'
ereturn  scalar eps =   `eps'
ereturn  scalar niter =  `k'
ereturn scalar widstat = e(widstat) 
ereturn scalar df_r = `dof'
ereturn scalar arf = e(arf)
ereturn local cmd "i2SLS_HDFE_MP"
ereturn local vcetype `option'
di in gr _col(55) "Number of obs = " in ye %8.0f e(N)
ereturn display	
	
}

********************************************************************************
*********************// ESTIMATION WITH FIXED EFFECTS //************************
********************************************************************************

if "`absorb'" != "" & "`ip'"==""{	
// case with covariates to absorb 
/*         CHECK FOR SEPERATION : CORREIRA CODE      */	
if "`nocheck'" == "1"{
loc tol = 1e-5
tempvar u w xb e
quietly: gen `u' =  !`depvar' if `touse'
quietly: su `u'  if `touse', mean
loc K = ceil(r(sum) / `tol' ^ 2)
quietly: gen `w' = cond(`depvar', `K', 1)  if `touse'
quietly: sum `w'
if r(mean)!=0{
while 1 {
quietly:	reghdfe `u' `_rhs'  `endog'  [fw=`w']  if `touse' , absorb(`absorb') resid(`e')
quietly:	predict double `xb'  if `touse', xbd
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
}
/*         DROP COLLINEAR VARIABLES      */	
	tempvar cste
	gen `cste' = 1
	quietly:   _rmcoll `_rhs' `endog' `cste' if `touse' , forcedrop 
	if r(k_omitted) >0 di 
	local alt_varlist `r(varlist)'
	local alt_varlist: list alt_varlist- endog
	local var_list `endog' `alt_varlist' 
	local instr_list `instr' `alt_varlist' 
	mata : delta = `delta'
/*         PREPARE iOLS      */	
	cap drop Z0_*
	cap drop E0_*
	cap drop M0_*
	cap drop Y0_*
	cap drop xb_hat*
	if "`alt_varlist'"=="" { // case with no X , only FE 
	quietly hdfe `endog' if `touse' [`weight'] , absorb(`absorb') generate(E0_)
	quietly hdfe `instr' if `touse' [`weight'] , absorb(`absorb') generate(Z0_)
	//quietly gen `y_tild' = log(max(1/`max_y',`depvar2')) if `touse'
	cap drop y_tild
quietly gen y_tild = ln(1+`depvar') if `touse'
	quietly	hdfe y_tild  if `touse' [`weight'] , absorb(`absorb') generate(Y0_)  tolerance(1e-3)  acceleration(sd) 
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
	quietly hdfe `alt_varlist'  if `touse'  [`weight'] , absorb(`absorb') generate(M0_)
	quietly hdfe `endog'  if `touse'  [`weight'] , absorb(`absorb') generate(E0_)
	quietly hdfe `instr'  if `touse'  [`weight'] , absorb(`absorb') generate(Z0_)
	cap drop y_tild  
	quietly gen y_tild = asinh(`depvar') if `touse'
	quietly	hdfe y_tild  if `touse'  [`weight'] , absorb(`absorb') generate(Y0_) tolerance(1e-3)  acceleration(sd)  
	local df_a = e(df_a)
	local dof_hdfe = e(df_a)
	mata : X=.
	mata : PX=.
	mata : PZ=.
	mata : y_tilde =.
	mata : Py_tilde =.
	mata : y =.
	mata : st_view(X,.,"`var_list'","`touse'")
	mata : st_view(PX,.,"E0_* M0_*","`touse'")
	mata : st_view(PZ,.,"Z0_* M0_*","`touse'")
	mata : st_view(y_tilde,.,"y_tild","`touse'")
	mata : st_view(Py_tilde,.,"Y0_","`touse'")
	mata : st_view(y,.,"`depvar'","`touse'")	
	}
/*         INITIAL VALUES      */	
mata : invPzX = invsym(cross(PX,PZ)*invsym(cross(PZ,PZ))*cross(PZ,PX))*cross(PX,PZ)*invsym(cross(PZ,PZ))
capture	 confirm matrix `from'
if _rc==0 {
	mata : beta_initial = st_matrix("`from'")
	mata : beta_initial = beta_initial'
}
else {
	mata : beta_initial = invPzX*cross(PZ,Py_tilde)
}
/*         PREPARE LOOP      */	
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
	local almost_conv = 1e-3
/*          LOOP      */	
mata: ivloop_function_fe("`touse'", y,xb_hat,xb_hat_M,PX,PZ,beta_initial,xb_hat_N,X,diff,Py_tilde,fe,y_tilde,delta,invPzX,beta_new,criteria,past_criteria)

/*         VARIANCE COVARIANCE CALCULATIONS      */	
	mata: ui = y:*exp(-xb_hat)
	mata: weight = ui:/(1 :+ delta)
	foreach var in `alt_varlist' {     // rename variables for last ols
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
quietly: ivreg2 Y0_ `alt_varlist' (`endog' = `instr') [`weight'`exp'] if `touse' , `option' noconstant   // standard case with X and FE 
local df_r = e(Fdf2) - `df_a'
 if "`cluster'" !="" {
 local df_r = e(Fdf2)
}
	foreach var in `alt_varlist' {      // rename variables back
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
/*         REPORT RESULTS       */	
	local names : colnames beta_final
	local nbvar : word count `names'
	mat rownames Sigma_tild = `names' 
    mat colnames Sigma_tild = `names' 
	local dof_final = e(df r)- `dof_hdfe'
	cap drop _COPY
	quietly: gen _COPY = `touse'
    ereturn post beta_final Sigma_tild , obs(`=e(N)') depname(`depvar') esample(`touse')  dof(`df_r')
	cap drop i2SLS_MP_HDFE_xb_hat
	cap drop i2SLS_MP_HDFE_fe
	cap drop i2SLS_MP_HDFE_error
	cap drop _reghdfe*
	cap drop y_tild
/*         RESTORE PRE-NORMALIZED VALUES       */	
//		mata: xb_hat = (xb_hat :+ ln(`max_y'))
		mata: st_store(., st_addvar("double", "i2SLS_MP_HDFE_fe"), "_COPY", fe)
	    mata: st_store(., st_addvar("double", "i2SLS_MP_HDFE_error"), "_COPY", ui)
    	mata: st_store(., st_addvar("double", "i2SLS_MP_HDFE_xb_hat"),"_COPY", xb_hat)
		cap drop _COPY
/*         EXPORT CONSTANTS       */	
ereturn scalar delta = `delta'
ereturn  scalar eps =   `eps'
ereturn  scalar niter =  `k'
ereturn local cmd "i2SLS_HDFE"
ereturn local vcetype `option'
di in gr _col(55) "Number of obs = " in ye %8.0f e(N)
ereturn display

* drop 
	cap drop E0_*
	cap drop Z0_*
	cap drop M0_* 
	cap drop Y0_*
	cap drop xb_hat*
}


********************************************************************************
**************// ESTIMATION WITH INCIDENTAL PARAMETER //************************
********************************************************************************


if "`absorb'" != "" & "`ip'"!=""{
// case with covariates to absorb 
/*         CHECK FOR SEPERATION : CORREIRA CODE      */	
if "`nocheck'" == "1"{
loc tol = 1e-5
tempvar u w xb e
quietly: gen `u' =  !`depvar' if `touse'
quietly: su `u'  if `touse', mean
loc K = ceil(r(sum) / `tol' ^ 2)
quietly: gen `w' = cond(`depvar', `K', 1)  if `touse'
quietly: sum `w'
if r(mean)!=0{
while 1 {
quietly:	reghdfe `u' `_rhs'  `endog'  [fw=`w']  if `touse' , absorb(`absorb') resid(`e')
quietly:	predict double `xb'  if `touse', xbd
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
}
/*         DROP COLLINEAR VARIABLES      */	
	tempvar cste
	gen `cste' = 1
	quietly:   _rmcoll `_rhs' `endog' `cste' if `touse' , forcedrop 
	if r(k_omitted) >0 di 
	local alt_varlist `r(varlist)'
	local alt_varlist: list alt_varlist- endog
	local var_list `endog' `alt_varlist' 
	local instr_list `instr' `alt_varlist' 
	mata : delta = `delta'
/*         PREPARE iOLS      */	
	cap drop Z0_*
	cap drop E0_*
	cap drop M0_*
	cap drop Y0_*
	cap drop xb_hat*
	if "`alt_varlist'"=="" { // case with no X , only FE 
	quietly hdfe `endog' if `touse' [`weight'] , absorb(`absorb') generate(E0_) acceleration(sd)   transform(sym)
	quietly hdfe `instr' if `touse' [`weight'] , absorb(`absorb') generate(Z0_) acceleration(sd)   transform(sym)
	//quietly gen `y_tild' = log(max(1/`max_y',`depvar2')) if `touse'
	cap drop y_tild
quietly gen y_tild = ln(1+`depvar') if `touse'
	quietly	hdfe y_tild  if `touse' [`weight'] , absorb(`absorb') generate(Y0_)  tolerance(1e-3)  acceleration(sd)   transform(sym)
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
	quietly hdfe `alt_varlist'  if `touse'  [`weight'] , absorb(`absorb') generate(M0_)
	quietly hdfe `endog'  if `touse'  [`weight'] , absorb(`absorb') generate(E0_)
	quietly hdfe `instr'  if `touse'  [`weight'] , absorb(`absorb') generate(Z0_)
	cap drop y_tild  
	quietly gen y_tild = ln(1+`depvar') if `touse'
	quietly	hdfe y_tild  if `touse'  [`weight'] , absorb(`absorb') generate(Y0_)  tolerance(1e-3)  acceleration(sd)   transform(sym) 
	local df_a = e(df_a)
	local dof_hdfe = e(df_a)
	mata : X=.
	mata : PX=.
	mata : PZ=.
	mata : y_tilde =.
	mata : Py_tilde =.
	mata : y =.
	mata : st_view(X,.,"`var_list'","`touse'")
	mata : st_view(PX,.,"E0_* M0_*","`touse'")
	mata : st_view(PZ,.,"Z0_* M0_*","`touse'")
	mata : st_view(y_tilde,.,"y_tild","`touse'")
	mata : st_view(Py_tilde,.,"Y0_","`touse'")
	mata : st_view(y,.,"`depvar'","`touse'")	
	}
/*         INITIAL VALUES      */	
mata : invPzX = invsym(cross(PX,PZ)*invsym(cross(PZ,PZ))*cross(PZ,PX))*cross(PX,PZ)*invsym(cross(PZ,PZ))
capture	 confirm matrix `from'
if _rc==0 {
	mata : beta_initial = st_matrix("`from'")
	mata : beta_initial = beta_initial'
}
else {
	mata : beta_initial = invPzX*cross(PZ,Py_tilde)
}
/*         PREPARE LOOP      */	
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
	local almost_conv = 1e-3
/*          LOOP      */	
mata: ivloop_function_ip("`touse'",y,xb_hat,PZ,beta_initial,X,diff,Py_tilde,y_tilde,delta,invPzX,beta_new,criteria,past_criteria)

/*         VARIANCE COVARIANCE CALCULATIONS      */	
	mata: ui = y:*exp(-xb_hat)
	mata: weight = ui:/(1 :+ delta)
	foreach var in `alt_varlist' {     // rename variables for last ols
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
quietly: ivreg2 Y0_ `alt_varlist' (`endog' = `instr') [`weight'`exp'] if `touse' , `option' noconstant   // standard case with X and FE 
local df_r = e(Fdf2) - `df_a'
 if "`cluster'" !="" {
 local df_r = e(Fdf2)
}
	foreach var in `alt_varlist' {      // rename variables back
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
/*         REPORT RESULTS       */	
	local names : colnames beta_final
	local nbvar : word count `names'
	mat rownames Sigma_tild = `names' 
    mat colnames Sigma_tild = `names' 
	local dof_final = e(df r)- `dof_hdfe'
	cap drop _COPY
	quietly: gen _COPY = `touse'
    ereturn post beta_final Sigma_tild , obs(`=e(N)') depname(`depvar') esample(`touse')  dof(`df_r')
	cap drop i2SLS_MP_HDFE_xb_hat
	cap drop i2SLS_MP_HDFE_fe
	cap drop i2SLS_MP_HDFE_error
	cap drop _reghdfe*
	cap drop y_tild
/*         RESTORE PRE-NORMALIZED VALUES       */	
//		mata: xb_hat = (xb_hat :+ ln(`max_y'))
	    mata: st_store(., st_addvar("double", "i2SLS_MP_HDFE_error"), "_COPY", ui)
    	mata: st_store(., st_addvar("double", "i2SLS_MP_HDFE_xb_hat"),"_COPY", xb_hat)
		cap drop _COPY
/*         EXPORT CONSTANTS       */	
ereturn scalar delta = `delta'
ereturn  scalar eps =   `eps'
ereturn  scalar niter =  `k'
ereturn local cmd "i2SLS_HDFE"
ereturn local vcetype `option'
di in gr _col(55) "Number of obs = " in ye %8.0f e(N)
ereturn display

* drop 
	cap drop E0_*
	cap drop Z0_*
	cap drop M0_* 
	cap drop Y0_*
	cap drop xb_hat*
}


end

cap: mata: mata drop ivloop_function_nofe()
cap: mata: mata drop ivloop_function_fe() 
cap: mata: mata drop ivloop_function_ip() 

mata:
void function ivloop_function_nofe(y,X,Z,beta_initial,delta,invPzX,criteria,xb_hat,y_tilde,beta_new,past_criteria)
{
max = strtoreal(st_local("maximum"))
lim = strtoreal(st_local("limit"))
show = (st_local("show"))
	for (i=1; i<=max;i++) {
	beta_initial[(cols(X)),1] = ln(mean(y:*exp(-X[.,1..(cols(X)-1)]*beta_initial[1..(cols(X)-1),1])))
	 xb_hat = X*beta_initial
	 y_tilde = ((y:*exp(-xb_hat) :- 1):/(1:+delta)) + xb_hat  
	 beta_new = invPzX*cross(Z,y_tilde)
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
void function ivloop_function_fe(string scalar touse, y,xb_hat,xb_hat_M,PX,PZ,beta_initial,xb_hat_N,X,diff,Py_tilde,fe,y_tilde,delta,invPzX,beta_new,criteria,past_criteria)
{
	max = strtoreal(st_local("maximum"))
lim = strtoreal(st_local("limit"))
show = (st_local("show"))
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
    stata("quietly: hdfe y_tild if \`touse' , absorb(\`absorb') generate(Y0_)  tolerance(\`almost_conv')  acceleration(sd)   transform(sym)")  
	st_view(Py_tilde,.,"Y0_",touse)
	beta_new = invPzX*cross(PZ,Py_tilde)
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
void function ivloop_function_ip(string scalar touse, y,xb_hat,PZ,beta_initial,X,diff,Py_tilde,y_tilde,delta,invPzX,beta_new,criteria,past_criteria)
{
max = strtoreal(st_local("maximum"))
lim = strtoreal(st_local("limit"))
show = (st_local("show"))
	for (i=1; i<=max;i++) {
	diff = y_tilde - Py_tilde
	xb_hat = X*beta_initial :+ ln(mean(exp(-X*beta_initial):*y))
	y_tilde = y:*exp(-xb_hat):/(1 :+ delta) + xb_hat  
	stata("cap drop y_tild")
	st_store(., st_addvar("double", "y_tild"), touse, y_tilde-diff)
	stata("cap drop Y0_")
    stata("quietly: hdfe y_tild if \`touse' , absorb(\`absorb') generate(Y0_)   tolerance(\`almost_conv')  acceleration(sd)   transform(sym) ")  
	st_view(Py_tilde,.,"Y0_",touse)
	beta_new = invPzX*cross(PZ,Py_tilde)
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

 use "C:\Users\ldpap\Downloadslog_gravity\Log of Gravity.dta" ,replace
  
  
