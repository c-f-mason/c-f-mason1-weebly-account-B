/*==============================================================================
  solar_jump.do
  Stata 19 translation of the GAUSS program "solar_jump_comand" (Solar data
  analysis, 02/22/17): constant-volatility models, NO GARCH.

  Model 1  (GAUSS proc lpr)       Geometric Brownian motion / plain normal
      y_t ~ N(mu, sigma^2)
  Model 2  (GAUSS proc mtlpjpr)   Mixed jump-diffusion, Poisson jump count
      y_t | K_t = k ~ N(mu + k*theta, sigma^2 + k*del^2),  K_t ~ Poisson(lambda)
      sum truncated at K = 10 jumps per period.

  NOTE on the GAUSS file as saved: it executes Model 1 and then hits `stop;`
  (line 102), so Model 2 is never run.  This translation has no `stop`: it runs both
  (set GJ_runjump = 0 to skip Model 2).

  Model 1 has a closed-form MLE (sample mean; sd with divisor N), so it is a
  built-in test of the machinery: the ml estimates must match it.

  Parameter transforms (Stata's ml has no bounds; GAUSS bounds in brackets):
      sigma = exp(.)          [0,100 / 0,1000]  upper bounds not imposed, checked after
      del   = exp(.)          [0,100]           upper bound not imposed, checked after
      lambda= invlogit(.)     [0,1]             open interval (0,1)
      mu, theta unrestricted  [-100,100]        bounds not imposed, checked after

  Untested against GAUSS output.  Each ml call sums the mixture in logs
  (log-sum-exp), mathematically identical to the GAUSS expression.
==============================================================================*/

version 19
clear all
macro drop GJ_*
set more off

*-------------------------------------------------------------------------------
* User settings
*-------------------------------------------------------------------------------
global GJ_datafile "/Users/chuckmason/Dropbox/Research/NeilWilmot_jumps/SREC/PJM_SREC_prices.dta"
global GJ_var      SRECp_ret        // variable to analyse
global GJ_scale    1                // oilp = GJ_scale * GJ_var.  1 = raw returns (matches the coauthor's numbers);
                                    // 100 = percent (the Stata preamble in the GAUSS file).  The MLE is scale-equivariant:
                                    // mu, sigma, theta, del scale by GJ_scale, lambda does not.
global GJ_maxopts  "difficult iterate(100) showtolerance"
global GJ_K        10               // maximum number of jumps per period
global GJ_multistart 1              // 1 = screen a grid of starts for Model 2 (GAUSS used one)
global GJ_runjump    1              // 1 = also estimate Model 2
global GJ_dropzero   0              // 1 = drop rows with an exact zero return (robustness check)
global GJ_droptails  0            // 1 = drop first AND last usable observations (after missing returns are removed)
                                  // 2 = drop first and last ROWS of the file (before missing returns are removed)
                                  // 3 = drop only the FIRST usable observation
                                  // 4 = drop only the LAST usable observation
global GJ_dropzero   0              // 1 = drop rows with an exact zero return (robustness check)
global GJ_droptails  0              // 1 = drop the first and last USABLE observations (after missing returns are removed)
                                    // 2 = drop the first and last ROWS of the file (before missing returns are removed)
global GJ_keepfirst  0              // n>0 = keep only the first n rows, as GAUSS's load solmat[n,k] does (0 = keep all)

* GAUSS start values.  Model 1: mu=5, sigma=14.  Model 2: Mu Sigma Lambda Theta Del
global GJ_gbm_mu0  = 5.0
global GJ_gbm_sig0 = 14.0
global GJ_mu0    = -0.01
global GJ_sig0   = 2.50
global GJ_lam0   = 0.10
global GJ_theta0 = 0.20
global GJ_del0   = 10.00
global GJ_sc = ${GJ_scale}/100      // the start values above are in the GAUSS file's x100 units; rescaled by this factor

*-------------------------------------------------------------------------------
* Data (loaded directly, as in the GAUSS file's Stata preamble)
*-------------------------------------------------------------------------------
display as text "Solar Jump paper estimation"
display as text " Single Jump processes"

use "${GJ_datafile}", clear
gen double oilp = ${GJ_scale}*${GJ_var}
if "${GJ_droptails}" == "2" {
    drop in 1
    drop in l
}
drop if oilp == .
if "${GJ_droptails}" == "1" {
    drop in 1
    drop in l
}
if "${GJ_droptails}" == "3" drop in 1
if "${GJ_droptails}" == "4" drop in l
if "${GJ_dropzero}" == "1" drop if oilp == 0
if ${GJ_keepfirst} > 0 keep in 1/${GJ_keepfirst}
display as text "Observations used: " _N

gen long t = _n
tsset t

summarize oilp, detail

* Exact zeros matter for Model 2: see the warning below.
quietly count if oilp == 0
local nzero = r(N)
display as text _n "Exact zero returns: `nzero' of " _N " (" %5.1f 100*`nzero'/_N "%)"
if `nzero' > 0 {
    display as text "Rows with an exact zero return (check: unchanged price, or a missing value coded as 0?):"
    list t ${GJ_var} if oilp == 0, noobs
    display as error "NOTE: in principle the Model 2 likelihood is unbounded as sigma -> 0 when any y equals mu exactly"
    display as error "(the k=0 component's density there ~ 1/sigma).  With only a few zeros the spike is narrow (log divergence);"
    display as error "watch the sigma column in the Model 2 screening table.  sigma collapsing toward 0 is this problem, not a solution."
}

*-------------------------------------------------------------------------------
* Model 1: geometric Brownian motion (GAUSS proc lpr)
*-------------------------------------------------------------------------------
capture program drop gj_gbm_lf
program define gj_gbm_lf
    version 19
    args lnf mu lnsig
    quietly replace `lnf' = -0.5*ln(2*_pi) - `lnsig' - 0.5*(($ML_y1 - `mu')/exp(`lnsig'))^2
end

display as text _n "Model 1: Geometric Brownian Motion"
matrix b1 = (${GJ_gbm_mu0}*${GJ_sc}, ln(${GJ_gbm_sig0}*${GJ_sc}))
ml model lf gj_gbm_lf (mu: oilp = ) /lnsigma, title("Geometric Brownian motion")
ml init b1, copy
ml maximize, difficult

estimates store gbm

* closed-form check (ML sd divides by N, not N-1)
quietly summarize oilp if e(sample)
local mhat = r(mean)
quietly generate double _dev2 = (oilp - `mhat')^2 if e(sample)
quietly summarize _dev2
local shat = sqrt(r(mean))
drop _dev2
matrix bb = e(b)
local mml = bb[1,1]
local sml = exp(bb[1,2])
display as text _n "N = " _N "  (compare these to the coauthor's Model 1 output digit by digit: they fingerprint the sample)"
display as text "Closed-form MLE:  mu = " %14.9f `mhat' "   sigma = " %14.9f `shat'
display as text   "ml estimate:      mu = " %14.9f `mml'  "   sigma = " %14.9f `sml'
if abs(`mml' - `mhat') > 1e-3*max(1, abs(`mhat')) | abs(`sml' - `shat') > 1e-3*`shat' {
    display as error "MISMATCH between ml and closed form: do not trust the machinery until this is resolved"
}
else display as text "ml matches the closed form."

capture noisily nlcom (Mu: _b[mu:_cons]) (Sigma: exp(_b[/lnsigma]))

if "${GJ_runjump}" != "1" exit

*-------------------------------------------------------------------------------
* Model 2: mixed jump-diffusion, Poisson jump count (GAUSS proc mtlpjpr)
*-------------------------------------------------------------------------------
capture program drop gj_jump_lf
program define gj_jump_lf
    version 19
    args lnf mu lnsig lamlgt th lndel
    tempvar sig lam del mx ssum
    quietly {
        generate double `sig' = exp(`lnsig')
        generate double `lam' = invlogit(`lamlgt')
        generate double `del' = exp(`lndel')
        * log of the k-th Poisson-weighted normal component
        forvalues k = 0/$GJ_K {
            tempvar a`k'
            generate double `a`k'' = `k'*ln(`lam') - lnfactorial(`k')                  ///
                - 0.5*ln(`sig'^2 + `k'*`del'^2)                                        ///
                - 0.5*($ML_y1 - `mu' - `k'*`th')^2/(`sig'^2 + `k'*`del'^2)
        }
        generate double `mx' = `a0'
        forvalues k = 1/$GJ_K {
            replace `mx' = max(`mx', `a`k'')
        }
        generate double `ssum' = 0
        forvalues k = 0/$GJ_K {
            replace `ssum' = `ssum' + exp(`a`k'' - `mx')
        }
        replace `lnf' = -`lam' - 0.5*ln(2*_pi) + `mx' + ln(`ssum')
    }
end

capture program drop gj_jfit
program define gj_jfit
    version 19
    args maxopts quiet
    ml model lf gj_jump_lf (mu: oilp = ) /lnsigma /lamlgt /theta /lndel, ///
        title("Mixed jump-diffusion (Poisson)")
    ml init b0, copy
    if "`quiet'" == "quiet" {
        capture ml maximize, `maxopts' nolog
        scalar gj_rc = _rc
    }
    else {
        capture noisily ml maximize, `maxopts'
        scalar gj_rc = _rc
    }
end

display as text _n "Model 2: (Multi) Mixed-jump diffusion process"

* start vector on the unconstrained scale: mu, ln(sigma), logit(lambda), theta, ln(del)
matrix b0 = (${GJ_mu0}*${GJ_sc}, ln(${GJ_sig0}*${GJ_sc}), logit(${GJ_lam0}), ${GJ_theta0}*${GJ_sc}, ln(${GJ_del0}*${GJ_sc}))

if "${GJ_multistart}" == "1" {
    * Mixture likelihoods are multimodal: screen a grid of starts, keep the best.
    scalar gj_best = -1e300
    foreach th in -20 -5 0.2 5 20 {
        foreach dl in 10 40 {
            foreach lm in 0.05 0.2 {
                matrix b0 = (${GJ_mu0}*${GJ_sc}, ln(${GJ_sig0}*${GJ_sc}), logit(`lm'), `th'*${GJ_sc}, ln(`dl'*${GJ_sc}))
                gj_jfit "difficult iterate(60)" quiet
                local ll = e(ll)
                matrix bb = e(b)
                display as text "start (x100 units) theta=`th' del=`dl' lambda=`lm':  ll = " %11.4f `ll'   ///
                    "  sigma=" %8.3f exp(bb[1,2]) "  lambda=" %6.4f invlogit(bb[1,3])       ///
                    "  theta=" %8.3f bb[1,4] "  del=" %8.3f exp(bb[1,5]) "  rc=" gj_rc
                if !missing(`ll') & `ll' > gj_best {
                    scalar gj_best = `ll'
                    matrix bbest = e(b)
                }
            }
        }
    }
    display as text _n "Best log likelihood from screening: " %11.4f gj_best
    if gj_best > -1e299 {
        matrix b0 = bbest
    }
    else {
        display as error "all screening fits failed; using the GAUSS start vector"
    }
}

gj_jfit "${GJ_maxopts}"
if gj_rc != 0 {
    display as error _n "WARNING: ml did not converge (rc = " gj_rc ").  Treat the results below with suspicion."
}
estimates store jump
display as text _n "N = " e(N) ";  log likelihood = " %11.4f e(ll) "  (in the units of the data as scaled)"
display as text "Log likelihood expressed in raw-return units (comparable across GJ_scale): " %11.4f e(ll) + e(N)*ln(${GJ_scale})

*-------------------------------------------------------------------------------
* Natural-scale estimates, bound check (GAUSS bounds), delta-method SEs
*-------------------------------------------------------------------------------
matrix bj = e(b)
local cn : colfullnames e(b)
display as text _n "coefficient names: `cn'"

local gmu  = bj[1,1]
local gsig = exp(bj[1,2])
local glam = invlogit(bj[1,3])
local gth  = bj[1,4]
local gdel = exp(bj[1,5])
display as text _n "Natural-scale estimates and bound check"
display as text "  mu     = " %12.5f `gmu'  cond(abs(`gmu'/${GJ_sc}) >= 100, "   <-- outside GAUSS bound [-100,100]", "")
display as text "  sigma  = " %12.5f `gsig' cond(`gsig'/${GJ_sc} > 100, "   <-- outside GAUSS bound [0,100]", cond(`gsig'/${GJ_sc} < 0.05, "   <-- collapsing toward 0: likely degenerate likelihood", ""))
display as text "  lambda = " %12.5f `glam' cond(`glam' < 1e-3, "   <-- at 0: theta and del are then not identified", cond(`glam' > 1-1e-3, "   <-- at upper bound 1", ""))
display as text "  theta  = " %12.5f `gth'  cond(abs(`gth'/${GJ_sc}) >= 100, "   <-- outside GAUSS bound [-100,100]", "")
display as text "  del    = " %12.5f `gdel' cond(`gdel'/${GJ_sc} > 100, "   <-- outside GAUSS bound [0,100]", "")

capture noisily nlcom (Mu: _b[mu:_cons]) (Sigma: exp(_b[/lnsigma]))       ///
    (Lambda: invlogit(_b[/lamlgt])) (Theta: _b[/theta]) (Del: exp(_b[/lndel]))
