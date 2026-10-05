/*==============================================================================
  solar_garch.do
  Stata 19 translation of the GAUSS program "solar_garch" (Solar data analysis,
  02/22/17; data switch updated May 2023): plain Gaussian GARCH(1,1), no jumps.

  Model (GAUSS proc garch):

      y_t = mu + e_t,   e_t ~ N(0, h_t)
      h_t = kappa + alpha*e_{t-1}^2 + beta*h_{t-1},      h_1 = sample mean of (y-mu)^2

  Parameter order and names match GAUSS: Mu Kappa Beta Alpha.
  Constraints in GAUSS: kappa >= .0001, beta >= 0, alpha >= 0, beta+alpha <= .99 (here the cap is the
  setting GG_cap, default .99, so its influence can be tested).

  Lessons carried over from solar_garch_jump.do / solar_jump.do:
   * ml hands the evaluator ONLY the estimation sample, so the sample is NOT
     restricted with "if t>1".  Row 1 seeds the recursion and its likelihood
     contribution is set to 0 inside Mata (GAUSS: __weight[1] = 0).
   * SCALE: estimates are scale-equivariant (mu by c, kappa by c^2, beta and
     alpha unchanged).  GJ-style settings below; bounds and start values written
     in GAUSS units are rescaled automatically.
   * SAMPLE switches (GG_droptails, GG_dropzero, GG_keepfirst) reproduce the
     ways a coauthor's sample can differ from yours.
   * The Mata likelihood is checked against a plain-Stata recomputation, and
     against Stata's own arch command as an independent implementation.

  Not carried over from the GAUSS file: the analytic gradient proc grd (used by
  CML as _cml_GradProc; ml here uses numerical derivatives), the commented
  IGARCH constraint block, and the commented jump-model experiments.

  Differences from GAUSS/CML you should know about
   * CML enforces the bounds and beta+alpha<=cap directly; ml cannot, so they are
     imposed by smooth reparameterisation (see gg_natural()).  If a bound binds
     (beta or alpha near 0, persistence near the cap) the transformed parameter runs
     to +/-infinity and the delta-method SEs are unreliable; a bound check is
     printed at the end.
   * GAUSS rescales the weights by N/(N-1); not replicated.  SEs differ by a
     factor of about sqrt((N-1)/N).
   * Untested against GAUSS output.
==============================================================================*/

version 19
clear all
macro drop GG_*
set more off

*-------------------------------------------------------------------------------
* User settings
*-------------------------------------------------------------------------------
global GG_datafile "/Users/chuckmason/Dropbox/Research/NeilWilmot_jumps/SREC/PJM_SREC_prices.dta"
global GG_var      SRECp_ret      // variable to analyse  <-- SET THESE THREE
global GG_scale    1              // oilp = GG_scale * GG_var.  GAUSS used 10*resid (solar=4), 100*rtn (2,3), 1*price (1)
global GG_boundscale 10           // data scale at which the GAUSS bound kappa >= .0001 and the start values were written
global GG_sortvar  ""             // optional: date variable to sort by; the recursion needs time order

* sample switches, applied in this order
global GG_droptails  0            // 1 = drop first AND last usable observations (after missing returns are removed)
                                  // 2 = drop first and last ROWS of the file (before missing returns are removed)
                                  // 3 = drop only the FIRST usable observation
                                  // 4 = drop only the LAST usable observation
global GG_dropzero   0            // 1 = drop rows with an exact zero value
global GG_keepfirst  0            // n>0 = keep only the first n rows, as GAUSS's load solmat[n,k] does

* Active-bound option.  If the bound report at the end says beta/alpha is at 0 or persistence is at the cap,
* ml cannot converge (the logistic transform only approaches a bound).  Pin the bound and estimate the rest:
global GG_fix ""                  // "" = none;  "persist" = beta+alpha fixed at the cap;  "beta0" = beta fixed at 0;  "alpha0" = alpha fixed at 0

* Cap on persistence beta+alpha.  GAUSS used .99.  The cap is a modelling choice: if the estimate sits on it, the
* result depends on it, so test .995, .999, 1 (integrated GARCH), or a large value such as 2 (effectively uncapped;
* beta+alpha >= 1 means no finite unconditional variance).
global GG_cap 0.99

global GG_maxopts "difficult iterate(100) showtolerance"
global GG_multistart 1            // 1 = fit from several (beta, alpha) starts and keep the best
global GG_benchmark  1            // 1 = also fit Stata's own arch command and show it for comparison
global GG_selfcheck  1            // 1 = recompute the log likelihood with a plain Stata loop and compare

* GAUSS start values: Mu = sample mean | Kappa = .5 | Beta = .1 | Alpha = .1   (kappa in GG_boundscale units)
global GG_kappa0 = 0.5
global GG_beta0  = 0.1
global GG_alpha0 = 0.1

*-------------------------------------------------------------------------------
* Data (loaded directly)
*-------------------------------------------------------------------------------
display as text "Solar GARCH(1,1) estimation (no jumps)"

use "${GG_datafile}", clear
if "${GG_sortvar}" != "" sort ${GG_sortvar}
gen double oilp = ${GG_scale}*${GG_var}
if "${GG_droptails}" == "2" {
    drop in 1
    drop in l
}
drop if oilp == .              // the recursion needs an unbroken series
if "${GG_droptails}" == "1" {
    drop in 1
    drop in l
}
if "${GG_droptails}" == "3" drop in 1
if "${GG_droptails}" == "4" drop in l
if "${GG_dropzero}" == "1" drop if oilp == 0
if ${GG_keepfirst} > 0 keep in 1/${GG_keepfirst}
display as text "Observations used (rows loaded): " _N

gen long t = _n                // row order is time order (check this!)
tsset t

summarize oilp, detail
local mbar = r(mean)

quietly count if oilp == 0
local nzero = r(N)
display as text _n "Exact zeros: `nzero' of " _N
if `nzero' > 0 {
    list t ${GG_var} if oilp == 0, noobs
}

global GG_Y oilp
global GG_sc = ${GG_scale}/${GG_boundscale}          // GAUSS-unit quantities are multiplied by this (variances by its square)
scalar gg_kmin = 0.0001*(${GG_sc})^2                  // GAUSS bound kappa >= .0001, in the data's units
if "${GG_cap}" == "" global GG_cap 0.99
scalar gg_cap = ${GG_cap}                             // cap on beta+alpha
global GG_mu0 = `mbar'                                // GAUSS start: Mu = meanc(y)

*-------------------------------------------------------------------------------
* Mata: parameter transformation and log-likelihood
*-------------------------------------------------------------------------------
* Unconstrained -> constrained:
*   kappa = kmin + exp(r2)                   (kappa >= kmin)
*   s = beta+alpha = cap*invlogit(r3)        (0 < s < cap)
*   w = beta's share = invlogit(r4)          beta = s*w, alpha = s*(1-w)
mata:
mata clear

// expand the free parameters to the full 4-vector when a bound is pinned (+/-40 = the bound to machine precision)
real rowvector gg_full(real rowvector r)
{
    string scalar fx
    fx = st_global("GG_fix")
    if (cols(r) == 4 | fx == "") return(r)
    if (fx == "persist") return( (r[1], r[2], 40, r[3]) )
    if (fx == "beta0")   return( (r[1], r[2], r[3], -40) )
    return( (r[1], r[2], r[3], 40) )                         // alpha0
}

real rowvector gg_natural(real rowvector r0)
{
    real rowvector r
    real scalar ub, kmin, s, w
    r    = gg_full(r0)
    ub   = st_numscalar("gg_cap")
    kmin = st_numscalar("gg_kmin")
    s    = ub*invlogit(r[3])
    w    = invlogit(r[4])
    return( (r[1], kmin + exp(r[2]), s*w, s*(1 - w)) )      // mu, kappa, beta, alpha
}

// observation-level log likelihood (GAUSS proc garch), stored in variable lnfvar
void gg_ll(string scalar lnfvar, string scalar yvar, real rowvector r)
{
    real rowvector b
    real colvector y, u2, h, lnf
    real scalar    t, n

    b = gg_natural(r)          // mu, kappa, beta, alpha

    y  = st_data(., yvar)      // ALL rows: ml must not subset the data
    n  = rows(y)
    u2 = (y :- b[1]):^2

    // GAUSS: recserar(kappa + lag(u2)*alpha, meanc(u2), beta)
    h    = J(n, 1, .)
    h[1] = mean(u2)
    for (t = 2; t <= n; t++) h[t] = b[2] + b[4]*u2[t-1] + b[3]*h[t-1]

    lnf = -0.5*(u2:/h :+ ln(2*pi()) :+ ln(h))

    // GAUSS: obs 1 has weight 0 (it seeds the recursion but is not in the likelihood)
    lnf[1] = 0

    st_store(., lnfvar, lnf)
}

// flag estimates sitting on (or near) a bound
void gg_bound_report(real rowvector r)
{
    real rowvector b
    real scalar    tol
    b   = gg_natural(r)
    tol = 0.001
    printf("\n{txt}Natural-scale estimates and bound check (data units as scaled)\n")
    printf("  mu     = %14.8f\n", b[1])
    printf("  kappa  = %14.8f%s\n", b[2], (b[2] < st_numscalar("gg_kmin")*(1 + tol) ? "   <-- at lower bound" : ""))
    printf("  beta   = %14.8f%s\n", b[3], (b[3] < tol ? "   <-- at 0 bound" : ""))
    printf("  alpha  = %14.8f%s\n", b[4], (b[4] < tol ? "   <-- at 0 bound" : ""))
    printf("  beta+alpha = %10.6f%s\n", b[3] + b[4], (b[3] + b[4] > st_numscalar("gg_cap") - tol ? "   <-- at the cap (" + strofreal(st_numscalar("gg_cap")) + ")" : ""))
    printf("  unconditional variance kappa/(1-beta-alpha) = %14.8f\n", b[2]/(1 - b[3] - b[4]))
    if (b[3] < tol | b[4] < tol | b[3] + b[4] > st_numscalar("gg_cap") - tol | b[2] < st_numscalar("gg_kmin")*(1 + tol)) {
        printf("{err}  A bound is active: the Hessian-based SEs are not valid.\n")
    }
}
end

*-------------------------------------------------------------------------------
* ml evaluator (method gf0: observation-level log likelihood, recursive in t)
*-------------------------------------------------------------------------------
capture program drop gg_eval
program define gg_eval
    version 19
    args todo b lnfj
    local nfree = 4
    if "$GG_fix" != "" local nfree = 3
    tempname p1 p2 p3 p4
    matrix gg_pv = J(1, `nfree', .)
    forvalues i = 1/`nfree' {
        mleval `p`i'' = `b', eq(`i') scalar
        matrix gg_pv[1, `i'] = `p`i''
    }
    mata: gg_ll("`lnfj'", "$GG_Y", st_matrix("gg_pv"))
end

*-------------------------------------------------------------------------------
* Guard against a partial run
*-------------------------------------------------------------------------------
if "${GG_Y}" == ""    global GG_Y oilp
if "${GG_mu0}" == ""  display as error "GG_* settings are not set: run the WHOLE do-file (do solar_garch.do), not a selection"
if "${GG_mu0}" == ""  exit 198

*-------------------------------------------------------------------------------
* Helpers: start vector on the unconstrained scale, and one ml fit
*-------------------------------------------------------------------------------
* kappa is given in GG_boundscale units; mu in the data's units.
capture program drop gg_startvec
program define gg_startvec
    version 19
    args mu kappa beta alpha
    local ub = scalar(gg_cap)
    local sum = `beta' + `alpha'
    local s0 = min(`sum', 0.9*`ub')       // keep the start inside the cap
    local sc = ${GG_sc}
    matrix b0 = ( `mu', ln(`kappa'*`sc'^2 - gg_kmin),     ///
                  logit(`s0'/`ub'), logit(`beta'/`sum') )
    if "${GG_fix}" == "persist" matrix b0 = (b0[1,1], b0[1,2], b0[1,4])
    if "${GG_fix}" == "beta0" | "${GG_fix}" == "alpha0" matrix b0 = (b0[1,1], b0[1,2], b0[1,3])
end

* One ml fit from matrix b0.  ml maximize aborts the do-file (r(430)) when it
* does not converge; capture so the script can carry on.  Sets scalar gg_rc.
capture program drop gg_fit
program define gg_fit
    version 19
    args maxopts quiet
    if "${GG_fix}" == "" {
        ml model gf0 gg_eval (mu: oilp = ) /lnkappa /apers /ashare, title("Gaussian GARCH(1,1)")
    }
    else if "${GG_fix}" == "persist" {
        ml model gf0 gg_eval (mu: oilp = ) /lnkappa /ashare, title("Gaussian GARCH(1,1), beta+alpha pinned at the cap")
    }
    else {
        local fxname "beta"
        if "${GG_fix}" == "alpha0" local fxname "alpha"
        ml model gf0 gg_eval (mu: oilp = ) /lnkappa /apers, title("Gaussian GARCH(1,1), `fxname' pinned at 0")
    }
    ml init b0, copy
    if "`quiet'" == "quiet" {
        capture ml maximize, `maxopts' nolog
        scalar gg_rc = _rc
    }
    else {
        capture noisily ml maximize, `maxopts'
        scalar gg_rc = _rc
    }
end

*-------------------------------------------------------------------------------
* Estimation.  GAUSS: CML, Newton-Raphson with step halving, analytic gradient.
* Here: ml's Newton-Raphson with numerical derivatives and 'difficult'.
* Row 1 seeds the recursion with zero weight (lnf[1] = 0 in Mata).
*-------------------------------------------------------------------------------
display as text _n "GARCH(1,1) estimation"

* GAUSS start vector: Mu = mean | Kappa | Beta | Alpha
gg_startvec ${GG_mu0} ${GG_kappa0} ${GG_beta0} ${GG_alpha0}

if "${GG_multistart}" == "1" {
    scalar gg_best = -1e300
    matrix bn = J(1, 4, .)
    foreach bt_al in "0.1 0.1" "0.05 0.3" "0.6 0.2" "0.3 0.6" {
        local bt : word 1 of `bt_al'
        local al : word 2 of `bt_al'
        gg_startvec ${GG_mu0} ${GG_kappa0} `bt' `al'
        gg_fit "difficult iterate(60)" quiet
        local ll = e(ll)
        capture mata: st_matrix("bn", gg_natural(st_matrix("e(b)")))
        display as text "start beta=`bt' alpha=`al':  ll = " %12.5f `ll'                    ///
            "  mu=" %10.6f bn[1,1] "  kappa=" %10.6f bn[1,2] "  beta=" %7.4f bn[1,3]      ///
            "  alpha=" %7.4f bn[1,4] "  rc=" gg_rc
        if !missing(`ll') & `ll' > gg_best {
            scalar gg_best = `ll'
            matrix bbest = e(b)
        }
    }
    display as text _n "Best log likelihood from screening: " %12.5f gg_best
    if gg_best > -1e299 {
        matrix b0 = bbest
    }
    else {
        display as error "all screening fits failed; using the GAUSS start vector"
    }
}

gg_fit "${GG_maxopts}"
if gg_rc != 0 {
    display as error _n "WARNING: ml did not converge (rc = " gg_rc ").  Treat the results below with suspicion;"
    display as error "the bound report at the end shows whether a constraint is active.  If it says beta or alpha is at 0, or"
    display as error "beta+alpha is at the cap, set GG_fix to beta0, alpha0 or persist at the top and rerun."
}

estimates store garch
scalar ll_garch = e(ll)
scalar n_garch  = e(N) - 1                 // the seed row contributes nothing
matrix braw = e(b)                         // raw (unconstrained) estimates, kept for the checks and bound report

display as text _n "N in the likelihood = " n_garch " (ml reports " (n_garch+1) " rows incl. the seed row);  log likelihood = " %12.5f ll_garch "  (data units as scaled)"
display as text "Log likelihood in raw units (comparable across GG_scale): " %12.5f (ll_garch + n_garch*ln(${GG_scale}))

*-------------------------------------------------------------------------------
* Self-check: recompute the log likelihood at the estimates with a plain Stata
* loop (no Mata) and compare with ml's value.
*-------------------------------------------------------------------------------
if "${GG_selfcheck}" == "1" {
    mata: st_matrix("bnat", gg_natural(st_matrix("braw")))
    local cmu  = bnat[1,1]
    local ckap = bnat[1,2]
    local cbet = bnat[1,3]
    local calp = bnat[1,4]

    quietly {
        generate double _chk_u2 = (oilp - `cmu')^2
        summarize _chk_u2
        generate double _chk_h = r(mean) in 1
        forvalues i = 2/`=_N' {
            replace _chk_h = `ckap' + `calp'*_chk_u2[`i'-1] + `cbet'*_chk_h[`i'-1] in `i'
        }
        generate double _chk_ll = -0.5*(_chk_u2/_chk_h + ln(2*_pi) + ln(_chk_h))
        summarize _chk_ll if t > 1
    }
    local llchk = r(sum)

    quietly {
        generate double _chk_lm = .
        mata: gg_ll("_chk_lm", "oilp", st_matrix("braw"))
        summarize _chk_lm if t > 1
    }
    local llm = r(sum)

    display as text _n "Self-check, three numbers that should agree:"
    display as text "  ml e(ll)                          = " %14.6f ll_garch
    display as text "  Mata gg_ll summed over t>1        = " %14.6f `llm'
    display as text "  plain-Stata loop summed over t>1  = " %14.6f `llchk'
    if abs(ll_garch - `llm') > 1e-6*max(1, abs(ll_garch)) | abs(`llm' - `llchk') > 1e-6*max(1, abs(`llm')) {
        display as error "SELF-CHECK FAILED: the three numbers disagree."
    }
    else display as text "Self-check passed."
    drop _chk_*
}

*-------------------------------------------------------------------------------
* Benchmark: Stata's own arch command, an independent implementation.  It
* initialises the variance differently from the GAUSS recursion and includes
* every row in the likelihood, so expect close but not identical estimates, and
* a different log likelihood.
*-------------------------------------------------------------------------------
if "${GG_benchmark}" == "1" {
    display as text _n "Benchmark: Stata arch GARCH(1,1) (all rows)"
    capture noisily arch oilp, arch(1) garch(1) nolog
    if _rc == 0 {
        matrix ba = e(b)
        display as text "arch coefficient vector (const, ARCH L1, GARCH L1, variance const):"
        matrix list ba, format(%12.6f) noheader
        display as text "(ours: mu, then below kappa/beta/alpha; arch's ARCH term = our alpha, its GARCH term = our beta, its variance constant = our kappa)"
    }
    estimates restore garch
}

*-------------------------------------------------------------------------------
* Results on the original (GAUSS) parameter scale, delta-method SEs
* (GAUSS: print b; print sqrt(diag(h)))
*-------------------------------------------------------------------------------
local km = gg_kmin
local cap = gg_cap
local S "(`cap'*invlogit(_b[/apers]))"
local W "invlogit(_b[/ashare])"
if "${GG_fix}" == "persist" local S "`cap'"
if "${GG_fix}" == "beta0"   local W "0"
if "${GG_fix}" == "alpha0"  local W "1"

local cn : colfullnames e(b)
display as text _n "coefficient names: `cn'"

capture noisily nlcom (Mu:      _b[mu:_cons])                           ///
      (Kappa:   `km' + exp(_b[/lnkappa]))                               ///
      (Beta:    `S'*`W')                                                ///
      (Alpha:   `S'*(1 - `W'))                                          ///
      (Persist: `S')                                                    ///
      (UncondVar: (`km' + exp(_b[/lnkappa]))/(1 - `S')), post

mata: gg_bound_report(st_matrix("braw"))

display as text _n "Units: data scaled by " ${GG_scale} ".  To compare with estimates in other units: mu scales by the ratio,"
display as text "kappa and the unconditional variance by its square; beta and alpha are unit-free."
