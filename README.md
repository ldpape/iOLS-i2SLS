# Iterated Ordinary Least Squares (iOLS) and Two Stage Least Squares (i2SLS) with High Dimensional Fixed Effects

This repository includes code for iOLS_MP_HDFE and i2SLS_MP_HDFE as described in Bellégo, Benatia and Pape (2021) : https://arxiv.org/abs/2203.11820

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
Without fixed effects, the package estimates a pseudo-GML model. 

        sysuse auto.dta, replace
        eststo: xi: iOLS_MP_HDFE price mpg i.foreign, robust

With fixed effects, the estimates do not suffer from the incidental parameter problem. 

        eststo: iOLS_MP_HDFE price mpg, absorb(foreign) robust
        esttab 

It is normal that both modes will not yield the same point estimates.

### i2SLS_MP_HDFE : Endogenous covariates
The syntax is different for i2SLS_MP_HDFE, as shown in the following example:

        use "http://fmwww.bc.edu/RePEc/bocode/i/ivp_bwt.dta", replace
        eststo: i2SLS_MP_HDFE bw parity white male, endog(cigspreg) instr(edfwhite edmwhite incwhite cigtax88)
        eststo: i2SLS_MP_HDFE bw parity white, endog(cigspreg) instr(edfwhite edmwhite incwhite cigtax88) absorb(male) robust
        esttab 



