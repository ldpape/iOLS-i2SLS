# Iterated Ordinary Least Squares (iOLS) and Two Stage Least Squares (i2SLS) with High Dimensional Fixed Effects

This repository includes code for iOLS_MP_HDFE and i2SLS_MP_HDFE as described in Bellégo, Benatia and Pape (2021) : https://arxiv.org/abs/2203.11820 .

This code is preliminary and subject to frequent updates. 

## Installation 
Before using this package, you must install the following:

          ssc install hdfe
          ssc install reghdfe
          ssc install moremata
          ssc install ivreg2
          ssc install ftools
          ssc install ppml
          ssc install ppmlhdfe 
          ssc install ranktest
          ssc install gtools

To install the beta version into Stata, run the following (requires at least Stata 14) : 

          cap ado uninstall iOLS_i2SLS
          net install iOLS_i2SLS, from("https://raw.githubusercontent.com/ldpape/iOLS-i2SLS/main/")
          
## Examples 
This package is compatible with the estout package.


### iOLS_MP_HDFE : Exogenous covariates
Without fixed effects, the package estimates a pseudo-GML model: use the prefix "xi" along with "i." if you include categorical variables.

        sysuse auto.dta, replace
        eststo: xi: iOLS_MP_HDFE price mpg i.foreign, robust

Use the "absorb" option to absorb fixed effects.  Both specifications will not yield the same point estimates.

        eststo: iOLS_MP_HDFE price mpg, absorb(foreign) robust

Include the option "warm" to first identify an approximate solution (using iOLS-ẟ) which is then rendered exact (using iOLS-ρ).

        eststo: iOLS_MP_HDFE price mpg, absorb(foreign) robust warm 

You can select the path to follow to approximate the final solution using "delta_path".

        eststo: iOLS_MP_HDFE price mpg, absorb(foreign) robust warm delta_path(1 10 100)

You can then display the output and export it using the "estout" package.

          esttab

### i2SLS_MP_HDFE : Endogenous covariates
The syntax is different for i2SLS_MP_HDFE, as shown in the following example:

        use "http://fmwww.bc.edu/RePEc/bocode/i/ivp_bwt.dta", replace
        eststo: i2SLS_MP_HDFE bw parity white male, endog(cigspreg) instr(edfwhite edmwhite incwhite cigtax88) cluster(male)
        eststo: i2SLS_MP_HDFE bw parity white, endog(cigspreg) instr(edfwhite edmwhite incwhite cigtax88) absorb(male) robust
        esttab 

The same options to find an approximate solution are available for i2SLS_MP_HDFE.

        use "http://fmwww.bc.edu/RePEc/bocode/i/ivp_bwt.dta", replace
        eststo: i2SLS_MP_HDFE bw parity white male, endog(cigspreg) instr(edfwhite edmwhite incwhite cigtax88) warm 
        eststo: i2SLS_MP_HDFE bw parity white, endog(cigspreg) instr(edfwhite edmwhite incwhite cigtax88) absorb(male) cluster(male) warm delta_path(1 10 100 1000)
        esttab 

##  Options :
- "limit( )" indicates the convergence threshold : it is set at 1e-3. Reset it by adding limit(1e-4), for example.
- "maximum( )" indicates the maximum number of iterations in the final step : it is set at 10,000. Reset it by adding maximum(500), for example.
- "rho( )" indicates the speed of convergence while finding the exact solution. It is set by default to one. Reset it by adding rho(10), for example.
- "show" displays each iteration's max absolute deviation compared to the previous step.
- "cluster( )" can be used to indicate on which variables standard errors should be clustered.

Weights are not allowed at this point in time.

## Citation : 
Please cite the following article : https://arxiv.org/abs/2203.11820

@misc{bbp2022,
      title={Dealing with Logs and Zeros in Regression Models}, 
      author={Christophe Bellégo and David Benatia and Louis Pape},
      year={2022},
      eprint={2203.11820},
      archivePrefix={arXiv},
      primaryClass={econ.EM},
      url={https://arxiv.org/abs/2203.11820}, 
}
