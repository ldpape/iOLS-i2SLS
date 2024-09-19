{smcl}
{* *! version 1.0 22march2021}{...}
{vieweralsosee "[R] poisson" "help poisson"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "reghdfe" "help reghdfe"}{...}
{vieweralsosee "ppml" "help ppml"}{...}
{vieweralsosee "ppmlhdfe" "help ppmlhdfe"}{...}
{viewerjumpto "Syntax" "i2SLS_MP_HDFE##syntax"}{...}
{viewerjumpto "Description" "i2SLS_MP_HDFE##description"}{...}
{viewerjumpto "Citation" "i2SLS_MP_HDFE##citation"}{...}
{viewerjumpto "Authors" "i2SLS_MP_HDFE##contact"}{...}
{viewerjumpto "Examples" "i2SLS_MP_HDFE##examples"}{...}
{viewerjumpto "Stored results" "i2SLS_MP_HDFE##results"}{...}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:i2SLS_MP_HDFE}} {hline 2} Iterated Two Stage Ordinary Least Squares (i2SLS) estimation of GPML with High-Dimensional Fixed Effects {p_end}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd}{cmd:Introduction:} This program implements Two Stage Ordinary Least Squares for Gamma Pseudo Maximum Likelihood (GPML), as described by {browse "https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3444996":Bellego, Benatia, and Pape (2021)}. {cmd:i2SLS_MP_HDFE} addresses the problem of the log of zero by iteratively running the "reghdfe" function. The rate of convergence can be increased using the {opt delta(#)} option, which defaults to 1.

{pstd}{cmd:Fixed-Effects:} The package applies a within-transformation to remove high-dimensional fixed effects, relying on the HDFE package developed by {browse "http://scorreia.com/research/hdfe.pdf":Sergio Correia (2017)}. The syntax is based on reghdfe.

{pstd}{cmd:Note:} This program checks for separation using the method proposed by {browse "https://arxiv.org/pdf/1903.01633.pdf":Correia, Guimarães, and Zylkin (2019)}, dropping problematic observations. It is important to check sample consistency.

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:i2SLS_MP_HDFE} {depvar} [{exog}] {ifin} {cmd:,} {opt delta(#)} {opt endog(varlist)} {opt instr(varlist)} {opt absorb(varlist)} [{help i2SLS_MP_HDFE##options:options}] {p_end}

{synoptset 22}{...}
{synopthdr:Variables}
{synoptline}
{synopt:{it:depvar}} Dependent variable{p_end}
{synopt:{it:exog}} List of exogenous explanatory variables (excludes instruments){p_end}
{synopt:{it:endog}} List of endogenous explanatory variables{p_end}
{synopt:{it:instr}} List of instruments (excludes exogenous variables){p_end}
{synopt:{it:fixed-effects}} List of categorical variables to be differenced-out{p_end}
{synoptline}
{p2colreset}{...}

{marker opt_summary}{...}
{synoptset 22 tabbed}{...}
{synopthdr:Options}
{synoptline}
{synopt:{opt delta(#)}} Set a strictly positive constant for delta (defaults to 1){p_end}
{synopt:{opt endog(varlist)}} Specifies endogenous variables{p_end}
{synopt:{opt instr(varlist)}} Specifies instruments{p_end}
{synopt:{opt absorb(varlist)}} Categorical variables treated as fixed effects{p_end}
{synopt:{opt vce(vcetype)}} Specifies the variance-covariance estimator: robust or clustered{p_end}
{synopt:{opt limit(#)}} Convergence criterion based on mean squared difference (defaults to 1e-3){p_end}
{synopt:{opt maximum(#)}} Maximum number of iterations (defaults to 10,000){p_end}
{synopt:{opt show}} Displays maximum absolute deviations between iterations{p_end}
{synopt:{opt ip}} Uses transformation immune to the incidental parameter problem{p_end}
{synoptline}

{marker Post-Estimation}{...}
{title:Post-Estimation}

{pstd} This program generates the following outcome variables:
{phang2}(i) Without fixed effects: {cmd:i2SLS_MP_HDFE_xb_hat} gives {it:X'b} in {it:Y = exp(X'b)U}, and {cmd:i2SLS_MP_HDFE_error} gives {it:U}.{p_end}
{phang2}(ii) With fixed effects: {cmd:i2SLS_MP_HDFE_xb_hat} gives {it:X'b + fixed-effects}, and {cmd:i2SLS_MP_HDFE_error} gives {it:U}.{p_end}
{phang2}(iii) With the {opt ip} option: {cmd:i2SLS_MP_HDFE_xb_hat} excludes fixed effects (i.e., {it:x1'b} in {it:Y=exp(x1'b + fe)U}), and {cmd:i2SLS_MP_HDFE_error} gives {it:U}.{p_end}

{marker authors}{...}
{title:Authors}

{pstd} Christophe Bellego, David Benatia, Louis Pape{p_end}
{pstd}CREST - ENSAE - HEC Montréal - Ecole Polytechnique{p_end}
{pstd}Contact: {browse "mailto:louis.pape@polytechnique.edu":louis.pape@polytechnique.edu}{p_end}

{marker citation}{...}
{title:Citation}

{pstd}
Bellégo Christophe, Benatia David, and Pape Louis-Daniel, Dealing with Logs and Zeros in Regression Models (2019).{break}
Série des Documents de Travail n° 2019-13. Available at SSRN: https://ssrn.com/abstract=3444996

{pstd}
BibTex citation: {break}
@misc{bellego_benatia_pape_2019,{break}
title={Dealing with Logs and Zeros in Regression Models},{break}
journal={SSRN},{break}
author={Bellégo, Christophe and Benatia, David and Pape, Louis-Daniel},{break}
year={2019},{break}
month={Sep}}

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. use "http://fmwww.bc.edu/RePEc/bocode/i/ivp_bwt.dta", replace}{p_end}
{phang2}{cmd:. i2SLS_MP_HDFE bw parity white male, endog(cigspreg) instr(edfwhite edmwhite incwhite cigtax88)}{p_end}
{phang2}{cmd:. i2SLS_MP_HDFE bw parity white, endog(cigspreg) instr(edfwhite edmwhite incwhite cigtax88) absorb(male) robust}{p_end}
{phang2}{cmd:. i2SLS_MP_HDFE bw parity white, endog(cigspreg) instr(edfwhite edmwhite incwhite cigtax88) absorb(male) robust ip}{p_end}
{phang2}{cmd:. i2SLS_MP_HDFE bw parity white, endog(cigspreg) instr(edfwhite edmwhite incwhite cigtax88) absorb(male) ip delta(3)}{p_end}
{phang2}{cmd:. i2SLS_MP_HDFE bw parity white, endog(cigspreg) instr(edfwhite edmwhite incwhite cigtax88) absorb(male) delta(3) limit(1e-5) maximum(100)}{p_end}

{marker results}{...}
{title:Stored Results}

{pstd}{cmd:i2SLS_MP_HDFE} stores the following in {cmd:e()}: {p_end}

{synoptset 24 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:e(N)}} Number of observations{p_end}
{synopt:{cmd:e(sample)}} Marks the sample used for estimation{p_end}
{synopt:{cmd:e(df_r)}} Degrees of freedom{p_end}
