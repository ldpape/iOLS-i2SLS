# iOLS_delta : This repository includes all available STATA programs  of iterated Ordinary Least Squares

This repository includes code for iOLS_delta, i2SLS_delta, iOLS_MP as well as iOLS_HDFE and i2SLS_HDFE which allows for high-dimensional fixed effects. These programs rely on other programs : you will need to have them installed.

>ssc install hdfe
>ssc install reghdfe
>ssc install moremata
>ssc install ivreg2

To install this code into Stata, run the following (requires at least Stata 14) : 

>cap ado uninstall iOLS_delta

>net install iOLS_delta, from("https://raw.githubusercontent.com/ldpape/iOLS_delta/master/")

