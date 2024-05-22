# iOLS_delta : This repository includes all available STATA programs  of iterated Ordinary Least Squares

This repository includes code for iOLS_MP_HDFE, i2SLS_MP_HDFE, iOLS_MP_HDFE_test, i2SLS_MP_HDFE, ppmlhdfe_test , popular_fix_test as described in Bellégo, Benatia and Pape (2021).

1). These programs rely on other programs : you will need to have them installed.

>ssc install hdfe

>ssc install reghdfe  // If an error appears, try downloading the more recent version from http://scorreia.com/software/reghdfe/install.html

>ssc install moremata

>ssc install ivreg2

>ssc install ftools

>ssc install ppml

>ssc install ppmlhdfe 

>ssc install ranktest

To install this code into Stata, run the following (requires at least Stata 14) : 

>cap ado uninstall iOLS_i2SLS

>net install iOLS_i2SLS, from("https://raw.githubusercontent.com/ldpape/iOLS_i2SLS/master/")

Please feel free to contact me to report a bug or ask a question. 

Note, this code is provided as is and may include potential errors.  It has been tested for Stata version 16 but has also worked on Stata 14. 

