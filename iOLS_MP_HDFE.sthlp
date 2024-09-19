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
{viewerjumpto "Stored results" "iOLS_MP_HDFE##results"}{...}

{title:Title}

{p2colset 5 18 20 2}{...}
{p2col :{cmd:iOLS_MP_HDFE} {hline 2}} Iterated Ordinary Least Squares (iOLS) estimation of GPML with High-Dimensional Fixed Effects {p_end}

{p2colreset}{...}

{pstd}{cmd:Introduction} This program implements iterated Ordinary Least Squares for Gamma Pseudo Maximum Likelihood (GPML), as described by {browse "https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3444996":Bellego, Benatia, and Pape (2021)}. {cmd: iOLS_MP_HDFE} is a solution to the problem of the log of zero.  This method relies on running the "reghdfe" function iteratively. The method depends on the rate of convergence, which you can increase with the option ", delta(number)" (set to 1 initially).
{pstd}{cmd:Fixed-Effects:} This package takes a within-transformation to difference out high-dimensional fixed effects. To do so, it relies on the HDFE package (used in reghdfe) developed by {browse "http://scorreia.com/research/hdfe.pdf": Sergio Correia (2017)}. The syntax follows from reghdfe.

{pstd}{cmd:Note:} This program automatically checks for the presence of seperation, which would preclude the existence of estimates, using the method proposed by {browse "https://arxiv.org/pdf/1903.01633.pdf": Correia, Guimarães, and Zylkin (2019) }. This results in dropping problematic observations. It is worth checking sample consistency.

{marker syntax}{...}
{title:Syntax}

{p 8 15 2} {cmd:iOLS_MP_HDFE}
{depvar} [{indepvars}]
{it:if}  {cmd:,} delta(#) absorb({it:fixed-effects}) [{help iOLS_MP_HDFE##options:options}] {p_end}

{synoptset 22}{...}
{synopthdr: variables}
{synoptline}
{synopt:{it:depvar}} Dependent variable{p_end}
{synopt:{it:indepvars}} List of exogenous explanatory variables {p_end}
{synopt:{it:fixed-effects}} List of categorical variables which are to be differenced-out {p_end}
{synoptline}
{p2colreset}{...}

{marker opt_summary}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt absorb}{cmd:(}{help iOLS_MP_HDFE##absorb:absorb}{cmd:)}}{it:absorb} This option allows you to include the categorical variables to treat as fixed effects. }{p_end}
{synopt:{opt vce}{cmd:(}{help iOLS_MP_HDFE##opt_vce:vcetype}{cmd:)}}{it:vcetype} May be classical if unspecified (assuming homoskedasticity), {opt r:obust}, or vce({opt cl:uster} varlist) (allowing two- and multi-way clustering) }{p_end}
{synopt:{opt delta}{cmd:(}{help iOLS_MP_HDFE##delta:delta}{cmd:)}}{it:delta} is any strictly positive constant. Set to 1 if unspecified.  }{p_end}
{synopt:{opt limit}{cmd:(}{help iOLS_MP_HDFE##limit:limit}{cmd:)}} Choose convergence criteria in terms of mean squared difference between two set of paramter estimates between two iterations. Set to 1e-3 if unspecified. }{p_end}
{synopt:{opt maximum}{cmd:(}{help iOLS_MP_HDFE##maximum:maximum}{cmd:)}} Maximum number of iterations. Set to 10,000 if unspecified. }{p_end}
{synopt:{opt show}{cmd:(}{help iOLS_MP_HDFE##show:show}{cmd:)}} This option shows the maximum absolute deviations between two iterations.  Convergence is assured when this number decreases regularly down to the convergence criteria  (1e-3) }{p_end}
{synopt:{opt ip}{cmd:(}{help iOLS_MP_HDFE##ip:ip}{cmd:)}} This option uses the transformation which is immune to the incidental parameter problem.  }{p_end}

{marker Post-Estimation}{...}
{title:Post-Estimation}

{pstd} This program generates outcome variables. (i) If you have no fixed-effects, "iOLS_MP_HDFE_xb_hat" is "X'b" in "Y = exp(X'b)U" and "iOLS_MP_HDFE_error" is equal to "U".  (ii) If you include fixed-effects with the "absorb()" option, "iOLS_MP_HDFE_fe" are the fixed effects, "iOLS_MP_HDFE_error" is the error "U", and "iOLS_MP_HDFE_xb_hat" is "X'b + fixed-effects".  (iii) If you use the "ip" option, "iOLS_MP_HDFE_error" is equal to "U" and  "iOLS_MP_HDFE_xb_hat" is equal to the linear index  excluding the fixed effects (i.e, equal to "x1'b" in  "Y=exp(x1'b + fe)U".

{marker authors}{...}
{title:Authors}

{pstd} Christophe Bellego, David Benatia, Louis Pape {break}
CREST - ENSAE - HEC Montréal - Ecole Polytechnique {break}
Contact: {browse "mailto:louis.pape@polytechnique.edu":louis.pape@polytechnique.edu} {p_end}

{marker citation}{...}
{title:Citation}

{pstd}
Bellégo Christophe, Benatia David, and Pape Louis-Daniel, Dealing with Logs and Zeros in Regression Models (2019).
Série des Documents de Travail n° 2019-13.
Available at SSRN: https://ssrn.com/abstract=3444996 

or in BibTex :

@misc{bellego_benatia_pape_2019, title={Dealing with logs and zeros in Regression Models}, journal={SSRN}, author={Bellégo, Christophe and Benatia, David and Pape, Louis-Daniel}, year={2019}, month={Sep}} 


{marker examples}{...}
{title:Examples}

{p_end}
{hline}
{phang2}{cmd:. sysuse auto.dta, replace }{p_end}
{phang2}{cmd:. xi: iOLS_MP_HDFE price mpg i.foreign, robust  }{p_end}
{phang2}{cmd:. iOLS_MP_HDFE price mpg , absorb(foreign) robust  }{p_end}
{phang2}{cmd:. iOLS_MP_HDFE price mpg , absorb(foreign) ip robust  }{p_end}
{phang2}{cmd:. iOLS_MP_HDFE price mpg , absorb(foreign) ip robust delta(3)  }{p_end}
{phang2}{cmd:. iOLS_MP_HDFE price mpg , absorb(foreign) ip robust delta(3) limit(1e-5) maximum(100) }{p_end}
{hline}


{marker results}{...}
{title:Stored results}

{pstd}
{cmd:iOLS_MP_HDFE} stores the following in {cmd:e()}:

{synoptset 24 tabbed}{...}
{syntab:Scalars}
{synopt:{cmd:e(N)}} number of observations{p_end}
{synopt:{cmd:e(sample)}} marks the sample used for estimation {p_end}
{synopt:{cmd:e(df_r)}} is the degrees of freedom {p_end}


