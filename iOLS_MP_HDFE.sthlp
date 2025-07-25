{smcl}
{* *! version 1.0 22march2021}{...}
{vieweralsosee "[R] poisson" "help poisson"}{...}
{vieweralsosee "" "--"}{...}
{vieweralsosee "reghdfe" "help reghdfe"}{...}
{vieweralsosee "ppml" "help ppml"}{...}
{vieweralsosee "ppmlhdfe" "help ppmlhdfe"}{...}
{viewerjumpto "Syntax" "iOLS_MP_HDFE##syntax"}{...}
{viewerjumpto "Description" "iOLS_MP_HDFE##description"}{...}
{viewerjumpto "Citation" "iOLS_MP_HDFE##citation"}{...}
{viewerjumpto "Authors" "iOLS_MP_HDFE##contact"}{...}
{viewerjumpto "Examples" "iOLS_MP_HDFE##examples"}{...}
{viewerjumpto "Description" "iOLS_MP_HDFE##Testing"}{...}
{viewerjumpto "Stored results" "iOLS_MP_HDFE##results"}{...}

{title:Title}

{p2colset 5 18 20 2}
{p2col :{cmd:iOLS_MP_HDFE} {hline 2}} Iterated Ordinary Least Squares (iOLS) estimation of GPML with High-Dimensional Fixed Effects {p_end}
{p2colreset}

{marker description}{...}
{title:Description}

{pstd} {cmd:iOLS_MP_HDFE} implements iterated Ordinary Least Squares for Gamma Pseudo Maximum Likelihood (GPML), as described by {browse "https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3444996":Bellego, Benatia, and Pape (2021)}. 
It addresses the problem of the log of zero by iteratively running the {cmd:reghdfe} function. Convergence is controlled with the option {cmd:rho(#). Estimation takes place in two-steps: 1) iOLS_delta provides an approximate solution to form a warm start (option warm) for any delta_path, 2) iOLS_MP_HDFE provides an exact solution.{p_end}

{pstd} The program applies a within-transformation to difference out high-dimensional fixed effects using the HDFE package from {browse "http://scorreia.com/research/hdfe.pdf":Sergio Correia (2017)}. The syntax is consistent with {cmd:reghdfe}.

{pstd} The program also checks for separation issues, dropping problematic observations based on the method by {browse "https://arxiv.org/pdf/1903.01633.pdf":Correia, Guimarães, and Zylkin (2019)}.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:iOLS_MP_HDFE} {depvar} [{indepvars}] {ifin}, {cmd:delta(#)} {cmd:absorb({it:fixed-effects})} [{help iOLS_MP_HDFE##options:options}] {p_end}

{synoptset 22}{...}
{synopthdr:Variables}
{synoptline}
{synopt:{it:depvar}} Dependent variable {p_end}
{synopt:{it:indepvars}} Explanatory variables {p_end}
{synopt:{it:fixed-effects}} Categorical variables to be differenced out {p_end}
{synoptline}

{marker opt_summary}{...}
{title:Options}

{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt rho(#)}} Set a strictly positive constant for initial rho (defaults to 1){p_end}
{synopt:{opt delta_path(#)}} Set a series of delta for warm startup such as delta_path(1 10 100){p_end}
{synopt:{opt absorb}{cmd:(}{help iOLS_MP_HDFE##absorb:absorb}{cmd:)}} Categorical variables to treat as fixed effects {p_end}
{synopt:{opt vce}{cmd:(}{help iOLS_MP_HDFE##opt_vce:vcetype}{cmd:)}} Specify variance-covariance estimator: {opt robust}, {opt cluster} for clustering {p_end}
{synopt:{opt limit}{cmd:(}{help iOLS_MP_HDFE##limit:limit}{cmd:)}} Convergence criteria (default: 1e-3) {p_end}
{synopt:{opt maximum}{cmd:(}{help iOLS_MP_HDFE##maximum:maximum}{cmd:)}} Maximum iterations (default: 10,000) {p_end}
{synopt:{opt warm}{cmd:(}{help iOLS_MP_HDFE##warm:warm}{cmd:)}} Warm startup using iOLS_delta{p_end}

{marker postestimation}{...}
{title:Post-Estimation}

{pstd} The program generates outcome variables:{p_end}
{phang2}1. Without fixed effects: {cmd:iOLS_MP_HDFE_xb_hat} = {cmd:X'b} and {cmd:iOLS_MP_HDFE_error} = U.{p_end}
{phang2}2. With fixed effects:  {cmd:iOLS_MP_HDFE_error} =  U =  {cmd:y*exp(-M_d X'b)} .{p_end}. {p_end}

{marker citation}{...}
{title:Citation}

{pstd}
Bellégo, Christophe, Benatia, David, and Pape, Louis-Daniel, Dealing with Logs and Zeros in Regression Models (2019). 
{browse "https://ssrn.com/abstract=3444996":SSRN Working Paper Series n° 2019-13}. {p_end}

{pstd}BibTeX:{p_end}
@misc{bellego_benatia_pape_2019, title={Dealing with logs and zeros in Regression Models}, journal={SSRN}, author={Bellégo, Christophe and Benatia, David and Pape, Louis-Daniel}, year={2019}, month={Sep}}

{marker examples}{...}
{title:Examples}

{phang2}{cmd:. sysuse auto.dta, replace}{p_end}
{phang2}{cmd:. xi: iOLS_MP_HDFE price mpg i.foreign,  robust}{p_end}
{phang2}{cmd:. xi: iOLS_MP_HDFE price mpg i.foreign, warm robust}{p_end}
{phang2}{cmd:. iOLS_MP_HDFE price mpg, absorb(foreign)  robust}{p_end}
{phang2}{cmd:. iOLS_MP_HDFE price mpg, absorb(foreign) warm  robust}{p_end}
{phang2}{cmd:. iOLS_MP_HDFE price mpg, absorb(foreign) warm delta_path(1 10 100) robust}{p_end}

{marker results}{...}
{title:Stored Results}

{pstd}{cmd:iOLS_MP_HDFE} stores the following in {cmd:e()}:{p_end}

{synoptset 24 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:e(N)}} Number of observations {p_end}
{synopt:{cmd:e(sample)}} Sample used for estimation {p_end}
{synopt:{cmd:e(df_r)}} Degrees of freedom {p_end}

{marker authors}{...}
{title:Authors}

{pstd}Christophe Bellego, David Benatia, Louis Pape {break}
CREST - ENSAE - HEC Montréal - Télécom Paris {break}
Contact: {browse "mailto:louis.pape@telecom-paris.fr":louis.pape@telecom-paris.fr}{p_end}
