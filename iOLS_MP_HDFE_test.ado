cap program drop iOLS_MP_HDFE_test
program define iOLS_MP_HDFE_test, eclass 
syntax varlist [if] [in] [aweight pweight fweight iweight] [, DELta(real 1) k(real 1)  LOGIT NONparametric LIMit(real 1e-3) from(name) MAXimum(real 10000) absorb(string)] 
	marksample touse
	local list_var `varlist'
	* get depvar and indepvar
	gettoken depvar list_var : list_var
	gettoken indepvar list_var : list_var, p("(")
 quietly:  iOLS_MP_HDFE `varlist' `excluded' if `touse' , delta(`delta') limit(`limit') from(`from') maximum(`maximum') absorb(`absorb')
  quietly: replace `touse' = e(sample)
	matrix beta_hat = e(b)
	matrix var_cov_beta_hat = e(V)
	* rhs of test
    tempvar lhs
	quietly: gen `lhs' = iOLS_MP_HDFE_error
tempvar dep_pos
quietly: gen `dep_pos' = `depvar'>0 if `touse'  
********************************************************************************
*                            PROBABILITY MODEL 	            	     	       *
********************************************************************************
if "`absorb'" != "" {
	if  "`logit'" =="logit" {
di in red "Using Logit Probability Model"
		local vlist1
		foreach item of varlist `absorb' {
		local vlist1 `vlist1' i.`item'
		}
		quietly:	xi: logit `dep_pos' `indepvar' `vlist1' if `touse'
		tempvar p_hat_temp
		quietly: predict `p_hat_temp' if `touse', pr 
		cap drop lambda_stat
		quietly: gen lambda_stat = 1/`p_hat_temp' if `touse'
		quietly: reg `lhs' lambda_stat if `dep_pos' & `touse', nocons 
		}
	else{
di in red "Using Linear Probability Model"
		quietly: reghdfe `dep_pos'  `indepvar'  if `touse' , absorb(`absorb') resid 
		tempvar p_hat_temp
		quietly: predict `p_hat_temp' if `touse', xbd
		quietly:		replace `p_hat_temp' = min(`p_hat_temp',1)
		quietly:        replace `p_hat_temp' = max(`p_hat_temp',0)
		cap drop lambda_stat
		quietly: gen lambda_stat = 1/`p_hat_temp' if `touse'
		quietly: reg `lhs' lambda_stat if `dep_pos' & `touse', nocons 
		}
}
else{
 if  "`nonparametric'" =="" {
di in red "Using Logit Probability Model"
quietly: logit `dep_pos' `indepvar' if `touse'
tempvar p_hat_temp
quietly:predict `p_hat_temp' if `touse', pr 
cap drop lambda_stat
quietly: gen lambda_stat = 1/`p_hat_temp' if `touse'
quietly: reg `lhs' lambda_stat if `dep_pos' & `touse', nocons 
	}
	else{
di in red "kNN Discrimination Probability Model"
tempvar p_hat_temp p_hat_neg
if `k'==1 {
quietly: sum `touse' if `touse'
local k = floor(sqrt(r(N))) 
}
 quietly: _rmcoll `indepvar' if `touse', forcedrop 
 quietly: discrim knn `r(varlist)'  if `touse' , k(`k') group(`dep_pos') notable ties(nearest)     priors(proportional) mahalanob
quietly: predict `p_hat_neg' `p_hat_temp'  if `touse', pr
*quietly: mrunning  `dep_pos'   `indepvar'  if `touse' , nograph predict(`p_hat_temp')
quietly: _pctile `p_hat_temp', p(5)
local w1=max(r(r1),1e-5)
quietly: _pctile `p_hat_temp', p(95)
local w2=min(r(r1),1) 
cap drop lambda_stat
quietly: gen lambda_stat =1/`p_hat_temp' if `touse'
quietly: reg `lhs' lambda_stat if `dep_pos' & `touse' & inrange(`p_hat_temp',`w1',`w2'), nocons       	
	} 
}
	matrix b = e(b)
	local lambda = _b[lambda_stat]	
	cap drop lambda_stat
******************************************************************************
*                   Return the information to STATA output		     		 *
******************************************************************************
di ""
di as result "Lambda Statistic = " in ye `lambda'
di "Interpretation: If the model is correct, lambda should be close to one. Reject lambda far from 1."
ereturn post b
ereturn scalar lambda = `lambda'
end

