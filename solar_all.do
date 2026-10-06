/*==============================================================================
  solar_all.do
  Stata 19.  One do-file, three analyses of the same sample, all translated
  from the GAUSS programs of Solar data analysis (Chuck, Neil, Luca):

    1. PLAIN JUMP-DIFFUSION (no GARCH), with geometric Brownian motion (GBM) as
       its benchmark.                          GAUSS: solar_jump_comand
    2. PLAIN GAUSSIAN GARCH(1,1), no jumps.    GAUSS: solar_garch
    3. GARCH(1,1) WITH POISSON JUMPS.          GAUSS: solar_garch_jump

  After every model there is a line marked

      >>>>>>>>>> ESTIMATES STORE (placeholder k of 7): ... <<<<<<<<<<
      estimates store NAME

  which stores the natural-scale estimates (GAUSS parameterisation, delta-method
  standard errors) under NAME, so the models can be put in one table.  Edit the
  NAME if you wish.  Raw (transformed-scale) ml results are stored as raw_*.
  A table template is at the end of the file.  e(ll_ours), e(N_ll), e(k_par)
  hold each model's log likelihood, observations in the likelihood, and number of
  free parameters.

  Models in the GARCH family (sections 2 and 3) use observation 1 only to seed
  the variance recursion (zero weight, as in GAUSS), so they use N-1 observations;
  the plain jump-diffusion uses all N unless JD_skipseed = 1 (below).

  Everything here has been run piecemeal in earlier versions (solar_jump.do,
  solar_garch.do, solar_garch_jump.do) and agreed with the coauthor's GAUSS
  output to three digits.  THIS COMBINED FILE HAS NOT YET BEEN RUN: run it once
  and compare each section with the earlier standalone results before use.

  Lessons built in (details in the standalone files):
   * ml hands a gf0 evaluator ONLY the estimation sample, so the sample is never
     restricted with "if t>1" in the GARCH models; the seed row's likelihood
     contribution is set to 0 inside Mata.
   * Estimates are scale-equivariant (mu, theta, del, sigma by c; kappa by c^2;
     beta, alpha, lambda unchanged).  Start values are rescaled automatically.
   * Mixture likelihoods are multimodal: start-value screening is on by default.
   * The Gaussian GARCH has persistence beta+alpha > 1 at its unconstrained
     optimum; the GAUSS cap of .99 is binding.  Section 2 therefore also fits the
     GARCH with the cap removed, which is the fair baseline for the jump model.
   * Stata's arch (Gaussian and Student-t) is run as an independent benchmark.

  Open modelling question: in section 3 the variance h_t is driven by the RAW
  squared residual (y-mu)^2, so a jump day feeds into next-period variance at
  full size.  A jump-adjusted recursion is the next thing to try.

  RUN THE WHOLE FILE (do solar_all.do), not a selection: settings, Mata and
  programs must all be defined in the same run.
==============================================================================*/

version 19
clear all
macro drop SC_* JD_* GG_* GJ_*        // start from clean settings
set more off

*===============================================================================
* COMMON SETTINGS  (data, sample, which sections to run)
*===============================================================================
global SC_datafile "/Users/chuckmason/Dropbox/Research/NeilWilmot_jumps/SREC/PJM_SREC_prices.dta"
global SC_var      SRECp_ret      // variable to analyse                <-- SET THESE
global SC_scale    1              // oilp = SC_scale * SC_var.  1 = raw returns (matches the coauthor's numbers), 100 = percent
global SC_sortvar  ""             // optional: date variable to sort by; the GARCH recursions need time order

* sample switches, applied in this order
global SC_droptails  0            // 1 = drop first AND last usable observations (after missing returns are removed)
                                  // 2 = drop first and last ROWS of the file (before missing returns are removed)
                                  // 3 = drop only the FIRST usable observation
                                  // 4 = drop only the LAST usable observation
global SC_dropzero   0            // 1 = drop rows with an exact zero return
global SC_keepfirst  0            // n>0 = keep only the first n rows, as GAUSS's load solmat[n,k] does (0 = keep all)

* which sections to run
global SC_run_jd 1                // 1 = plain jump-diffusion (+ GBM)
global SC_run_gg 1                // 1 = plain Gaussian GARCH(1,1)
global SC_run_gj 1                // 1 = GARCH(1,1) with Poisson jumps

* shared estimation settings
global SC_maxopts "difficult iterate(100) showtolerance"   // options for ml maximize (defaults: tolerance(1e-6) ltolerance(1e-7) nrtolerance(1e-5))
global SC_multistart 1            // 1 = screen a grid of start values first, then polish the best (GAUSS used one start)
global SC_selfcheck  1            // 1 = recompute each GARCH log likelihood with a plain Stata loop and compare
global SC_K          10           // maximum number of jumps per period in the Poisson sums

*===============================================================================
* MODEL-SPECIFIC SETTINGS  (the defaults reproduce the GAUSS programs)
*===============================================================================
* ---- Section 1: plain jump-diffusion and GBM  (GAUSS start values, written in x100 units)
global JD_scale      ${SC_scale}
global JD_maxopts    "${SC_maxopts}"
global JD_K          ${SC_K}
global JD_multistart ${SC_multistart}
global JD_skipseed   0            // 1 = estimate on t>1 only, so the sample matches the GARCH models exactly
                                  // (0 keeps all rows, which reproduced the coauthor's jump-only results)
global JD_gbm_mu0  = 5.0
global JD_gbm_sig0 = 14.0
global JD_mu0    = -0.01
global JD_sig0   = 2.50
global JD_lam0   = 0.10
global JD_theta0 = 0.20
global JD_del0   = 10.00
global JD_sc = ${JD_scale}/100    // start values above are in x100 units
global JD_if "t > 0"
if "${JD_skipseed}" == "1" global JD_if "t > 1"

* ---- Section 2: plain Gaussian GARCH
global GG_scale      ${SC_scale}
global GG_boundscale 10           // data scale at which the GAUSS bound kappa >= .0001 and the start values were written
global GG_cap        0.99         // cap on persistence beta+alpha (GAUSS used .99; the cap BINDS on this data)
global GG_fix        ""           // "" = none; "persist" = beta+alpha fixed at the cap; "beta0" = beta fixed at 0; "alpha0" = alpha fixed at 0
                                  // (use when the bound report says a bound is active and ml will not converge)
global GG_extra_uncapped 1        // 1 = also fit the GARCH with the cap removed (the fair baseline for the jump model)
global GG_cap_u      5            // the "removed" cap (beta+alpha < 5 is effectively unrestricted here)
global GG_maxopts    "${SC_maxopts}"
global GG_multistart ${SC_multistart}
global GG_benchmark  1            // 1 = also fit Stata's own arch (Gaussian) and show it for comparison
global GG_benchmark_t 1           // 1 = also fit arch with Student-t innovations: heavy tails WITHOUT jumps
global GG_selfcheck  ${SC_selfcheck}
global GG_kappa0 = 0.5            // GAUSS start values: Mu = sample mean | Kappa | Beta | Alpha   (kappa in GG_boundscale units)
global GG_beta0  = 0.1
global GG_alpha0 = 0.1

* ---- Section 3: GARCH(1,1) with Poisson jumps
global GJ_scale      ${SC_scale}
global GJ_startscale 100          // data units in which the start values / screening grid below are written
global GJ_maxopts    "${SC_maxopts}"
global GJ_multistart ${SC_multistart}
global GJ_selfcheck  ${SC_selfcheck}
global GJ_K          ${SC_K}
global GJ_mu0 = 0.10              // GAUSS start values: Mu | Kappa | Beta | Alpha | Lambda | Theta | Del
global GJ_kappa0 = 1.5
global GJ_beta0 = 0.05
global GJ_alpha0 = 0.3
global GJ_lam0 = 0.3
global GJ_theta0 = 0.50
global GJ_del0 = 3.5

*===============================================================================
* DATA  (loaded once; all three analyses use exactly this sample)
*===============================================================================
display as text "Solar jump paper: plain jump-diffusion, plain GARCH, GARCH with jumps"

use "${SC_datafile}", clear
if "${SC_sortvar}" != "" sort ${SC_sortvar}
gen double oilp = ${SC_scale}*${SC_var}
if "${SC_droptails}" == "2" {
    drop in 1
    drop in l
}
drop if oilp == .              // the recursions need an unbroken series
if "${SC_droptails}" == "1" {
    drop in 1
    drop in l
}
if "${SC_droptails}" == "3" drop in 1
if "${SC_droptails}" == "4" drop in l
if "${SC_dropzero}" == "1" drop if oilp == 0
if ${SC_keepfirst} > 0 keep in 1/${SC_keepfirst}
display as text "Observations used (rows loaded): " _N

gen long t = _n                // row order is time order (check this!)
tsset t

summarize oilp, detail
local mbar = r(mean)

quietly count if oilp == 0
local nzero = r(N)
display as text _n "Exact zero returns: `nzero' of " _N
if `nzero' > 0 {
    display as text "Rows with an exact zero return (unchanged price, or a missing value coded as 0?):"
    list t ${SC_var} if oilp == 0, noobs
    display as error "NOTE: in principle the plain jump-diffusion likelihood is unbounded as sigma -> 0 when any y equals mu exactly."
    display as error "With only a few zeros the spike is narrow; watch the sigma column in the section-1 screening table."
}

* derived quantities used by the sections
global GG_Y oilp
global GJ_Y oilp
global GG_sc = ${GG_scale}/${GG_boundscale}           // GAUSS-unit quantities are multiplied by this (variances by its square)
scalar gg_kmin = 0.0001*(${GG_sc})^2                  // GAUSS bound kappa >= .0001, in the data's units
scalar gg_cap  = ${GG_cap}                            // cap on beta+alpha
global GG_mu0 = `mbar'                                // GAUSS start: Mu = meanc(y)
global GJ_sc  = ${GJ_scale}/${GJ_startscale}          // start values are multiplied by this (kappa by its square)

* results scalars, initialised to missing so the comparison table works even if a section is switched off
foreach m in gbm jd garch garch_u jump arch archt {
    scalar ll_`m' = .
    scalar n_`m'  = .
    scalar k_`m'  = .
}

*===============================================================================
* MATA: parameter transformations and log likelihoods for the two GARCH models
*===============================================================================
mata:
mata clear

// ---------------------------------------------------------------- plain GARCH (gg_*)

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

// Nelson (1990): GARCH(1,1) with Gaussian z is strictly stationary iff E ln(beta + alpha z^2) < 0
real scalar gg_lyap(real scalar b, real scalar a)
{
    real colvector z
    real scalar    n, L, h
    L = 12
    n = 24001
    h = 2*L/(n - 1)
    z = rangen(-L, L, n)
    return( h*sum( normalden(z) :* ln(max((b, 1e-12)) :+ a*z:^2) ) )
}

// same exponent when z is a Student-t with nu df, scaled to unit variance (nu > 2)
real scalar gg_lyap_t(real scalar b, real scalar a, real scalar nu)
{
    real colvector z
    real scalar    n, L, h, c
    if (missing(nu) | nu <= 2) return(.)
    L = 300
    n = 300001
    h = 2*L/(n - 1)
    c = sqrt(nu/(nu - 2))
    z = rangen(-L, L, n)
    return( h*sum( c*tden(nu, c*z) :* ln(max((b, 1e-12)) :+ a*z:^2) ) )
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
    if (b[3] + b[4] < 1) printf("  unconditional variance kappa/(1-beta-alpha) = %14.8f\n", b[2]/(1 - b[3] - b[4]))
    else printf("  unconditional variance: undefined (beta+alpha >= 1)\n")
    printf("  Nelson (1990) exponent E ln(beta+alpha z^2) = %9.5f   (> 0: not even strictly stationary)\n", gg_lyap(b[3], b[4]))
    if (b[3] < tol | b[4] < tol | b[3] + b[4] > st_numscalar("gg_cap") - tol | b[2] < st_numscalar("gg_kmin")*(1 + tol)) {
        printf("{err}  A bound is active: the Hessian-based SEs are not valid.\n")
    }
}

// ---------------------------------------------------------------- GARCH + jumps (gj_*)

real rowvector gj_natural(real rowvector r)
{
    real scalar lb, ub, s, w
    lb = 0.0001
    ub = 0.99
    s  = 2*lb + (ub - 2*lb)*invlogit(r[3])
    w  = invlogit(r[4])
    return( (r[1], exp(r[2]), lb + (s - 2*lb)*w, lb + (s - 2*lb)*(1 - w),
             lb + (1 - lb)*invlogit(r[5]), r[6], lb + exp(r[7])) )
}

// observation-level log likelihood (GAUSS proc garchmj), stored in variable lnfvar
void gj_ll(string scalar lnfvar, string scalar yvar, string rowvector pn,
           real scalar K)
{
    real rowvector r, b
    real colvector y, u2, h, v, m, lnf
    real matrix    A
    real scalar    i, k, t, n

    r = J(1, 7, .)
    for (i = 1; i <= 7; i++) r[i] = st_numscalar(pn[i])
    b = gj_natural(r)          // mu, kappa, beta, alpha, lambda, theta, del

    y  = st_data(., yvar)
    n  = rows(y)
    u2 = (y :- b[1]):^2

    // GAUSS: recserar(kappa + lag(u2)*alpha, meanc(u2), beta)
    h    = J(n, 1, .)
    h[1] = mean(u2)
    for (t = 2; t <= n; t++) h[t] = b[2] + b[4]*u2[t-1] + b[3]*h[t-1]

    // log of each Poisson-weighted normal component, k = 0..K
    A = J(n, K + 1, .)
    for (k = 0; k <= K; k++) {
        v = h :+ k*b[7]^2
        A[., k+1] = k*ln(b[5]) - lnfactorial(k) :- 0.5*ln(v) :-
                    0.5*((y :- b[1] :- k*b[6]):^2):/v
    }
    m   = rowmax(A)
    lnf = -b[5] - 0.5*ln(2*pi()) :+ m :+ ln(rowsum(exp(A :- m)))

    // GAUSS: obs 1 has weight 0 (it seeds the recursion but is not in the likelihood).
    // Done here, not with an "if t>1" sample: ml hands the evaluator ONLY the
    // estimation sample, which would silently start the recursion at row 2.
    lnf[1] = 0

    st_store(., lnfvar, lnf)
}

// flag estimates sitting on (or near) a GAUSS bound
void gj_bound_report(real rowvector r)
{
    real rowvector b
    real scalar    tol
    b   = gj_natural(r)
    tol = 0.001
    printf("\n{txt}Natural-scale estimates and bound check (data units as scaled)\n")
    printf("  mu     = %12.6f\n", b[1])
    printf("  kappa  = %12.6f\n", b[2])
    printf("  beta   = %12.6f%s\n", b[3], (b[3] < 0.0001 + tol ? "   <-- at lower bound" : ""))
    printf("  alpha  = %12.6f%s\n", b[4], (b[4] < 0.0001 + tol ? "   <-- at lower bound" : ""))
    printf("  lambda = %12.6f%s\n", b[5], (b[5] < 0.0001 + tol ? "   <-- at lower bound" :
                                          (b[5] > 1 - tol ? "   <-- at upper bound" : "")))
    printf("  theta  = %12.6f\n", b[6])
    printf("  del    = %12.6f%s\n", b[7], (b[7] < 0.0001 + tol ? "   <-- at lower bound" : ""))
    printf("  beta+alpha = %9.6f%s\n", b[3] + b[4], (b[3] + b[4] > 0.99 - tol ? "   <-- at .99 cap" : ""))
    if (b[3] < 0.0001 + tol | b[4] < 0.0001 + tol | b[3] + b[4] > 0.99 - tol |
        b[5] < 0.0001 + tol | b[5] > 1 - tol | b[7] < 0.0001 + tol) {
        printf("{err}  A bound is active: the Hessian-based SEs are not valid.\n")
    }
}
end

*===============================================================================
* PROGRAMS (all defined here, at the top level, before any model is run)
*===============================================================================
capture program drop jd_gbm_lf
program define jd_gbm_lf
    version 19
    args lnf mu lnsig
    quietly replace `lnf' = -0.5*ln(2*_pi) - `lnsig' - 0.5*(($ML_y1 - `mu')/exp(`lnsig'))^2
end

capture program drop jd_jump_lf
program define jd_jump_lf
    version 19
    args lnf mu lnsig lamlgt th lndel
    tempvar sig lam del mx ssum
    quietly {
        generate double `sig' = exp(`lnsig')
        generate double `lam' = invlogit(`lamlgt')
        generate double `del' = exp(`lndel')
        * log of the k-th Poisson-weighted normal component
        forvalues k = 0/$JD_K {
            tempvar a`k'
            generate double `a`k'' = `k'*ln(`lam') - lnfactorial(`k')                  ///
                - 0.5*ln(`sig'^2 + `k'*`del'^2)                                        ///
                - 0.5*($ML_y1 - `mu' - `k'*`th')^2/(`sig'^2 + `k'*`del'^2)
        }
        generate double `mx' = `a0'
        forvalues k = 1/$JD_K {
            replace `mx' = max(`mx', `a`k'')
        }
        generate double `ssum' = 0
        forvalues k = 0/$JD_K {
            replace `ssum' = `ssum' + exp(`a`k'' - `mx')
        }
        replace `lnf' = -`lam' - 0.5*ln(2*_pi) + `mx' + ln(`ssum')
    }
end

capture program drop jd_jfit
program define jd_jfit
    version 19
    args maxopts quiet
    ml model lf jd_jump_lf (mu: oilp = ) /lnsigma /lamlgt /theta /lndel if ${JD_if}, ///
        title("Mixed jump-diffusion (Poisson)")
    ml init b0, copy
    if "`quiet'" == "quiet" {
        capture ml maximize, `maxopts' nolog
        scalar jd_rc = _rc
    }
    else {
        capture noisily ml maximize, `maxopts'
        scalar jd_rc = _rc
    }
end

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

capture program drop gj_eval
program define gj_eval
    version 19
    args todo b lnfj
    tempname p1 p2 p3 p4 p5 p6 p7
    forvalues i = 1/7 {
        mleval `p`i'' = `b', eq(`i') scalar
    }
    mata: gj_ll("`lnfj'", "$GJ_Y", ("`p1'","`p2'","`p3'","`p4'","`p5'","`p6'","`p7'"), $GJ_K)
end

capture program drop gj_startvec
program define gj_startvec
    version 19
    args mu kappa beta alpha lam theta del
    local lb = 0.0001
    local ub = 0.99
    local s0 = `beta' + `alpha'
    local sc = ${GJ_sc}
    matrix b0 = ( `mu'*`sc', ln(`kappa'*`sc'^2),                  ///
                  logit((`s0' - 2*`lb')/(`ub' - 2*`lb')),         ///
                  logit((`beta' - `lb')/(`s0' - 2*`lb')),         ///
                  logit((`lam' - `lb')/(1 - `lb')),               ///
                  `theta'*`sc', ln(`del'*`sc' - `lb') )
end

capture program drop gj_fit
program define gj_fit
    version 19
    args maxopts quiet
    ml model gf0 gj_eval (mu: oilp = ) /lnkappa /apers /ashare /lamlgt /theta /lndel, ///
        title("GARCH(1,1) with Poisson jumps")
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


*===============================================================================
* SECTION 1.  PLAIN JUMP-DIFFUSION (no GARCH), with geometric Brownian motion as
*             its benchmark.  GAUSS: solar_jump_comand (procs lpr and mtlpjpr).
*===============================================================================
if "${SC_run_jd}" == "1" {
display as text _n "Model 1: Geometric Brownian Motion"
matrix b1 = (${JD_gbm_mu0}*${JD_sc}, ln(${JD_gbm_sig0}*${JD_sc}))
ml model lf jd_gbm_lf (mu: oilp = ) /lnsigma if ${JD_if}, title("Geometric Brownian motion")
ml init b1, copy
ml maximize, difficult

estimates store raw_gbm
scalar ll_gbm = e(ll)
scalar n_gbm  = e(N)
scalar k_gbm  = 2

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

* natural-scale estimates (delta-method SEs), ready for a table
estimates restore raw_gbm
capture noisily nlcom (Mu: _b[mu:_cons]) (Sigma: exp(_b[/lnsigma])), post
ereturn scalar ll_ours = ll_gbm
ereturn scalar N_ll    = n_gbm
ereturn scalar k_par   = k_gbm
* >>>>>>>>>> ESTIMATES STORE (placeholder 1 of 7): geometric Brownian motion <<<<<<<<<<
estimates store gbm


display as text _n "Model 2: (Multi) Mixed-jump diffusion process"

* start vector on the unconstrained scale: mu, ln(sigma), logit(lambda), theta, ln(del)
matrix b0 = (${JD_mu0}*${JD_sc}, ln(${JD_sig0}*${JD_sc}), logit(${JD_lam0}), ${JD_theta0}*${JD_sc}, ln(${JD_del0}*${JD_sc}))

if "${JD_multistart}" == "1" {
    * Mixture likelihoods are multimodal: screen a grid of starts, keep the best.
    scalar jd_best = -1e300
    foreach th in -20 -5 0.2 5 20 {
        foreach dl in 10 40 {
            foreach lm in 0.05 0.2 {
                matrix b0 = (${JD_mu0}*${JD_sc}, ln(${JD_sig0}*${JD_sc}), logit(`lm'), `th'*${JD_sc}, ln(`dl'*${JD_sc}))
                jd_jfit "difficult iterate(60)" quiet
                local ll = e(ll)
                matrix bb = e(b)
                display as text "start (x100 units) theta=`th' del=`dl' lambda=`lm':  ll = " %11.4f `ll'   ///
                    "  sigma=" %8.3f exp(bb[1,2]) "  lambda=" %6.4f invlogit(bb[1,3])       ///
                    "  theta=" %8.3f bb[1,4] "  del=" %8.3f exp(bb[1,5]) "  rc=" jd_rc
                if !missing(`ll') & `ll' > jd_best {
                    scalar jd_best = `ll'
                    matrix bbest = e(b)
                }
            }
        }
    }
    display as text _n "Best log likelihood from screening: " %11.4f jd_best
    if jd_best > -1e299 {
        matrix b0 = bbest
    }
    else {
        display as error "all screening fits failed; using the GAUSS start vector"
    }
}

jd_jfit "${JD_maxopts}"
if jd_rc != 0 {
    display as error _n "WARNING: ml did not converge (rc = " jd_rc ").  Treat the results below with suspicion."
}
estimates store raw_jd
scalar ll_jd = e(ll)
scalar n_jd  = e(N)
scalar k_jd  = 5
display as text _n "N = " e(N) ";  log likelihood = " %11.4f e(ll) "  (in the units of the data as scaled)"
display as text "Log likelihood expressed in raw-return units (comparable across JD_scale): " %11.4f e(ll) + e(N)*ln(${JD_scale})

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
display as text "  mu     = " %12.5f `gmu'  cond(abs(`gmu'/${JD_sc}) >= 100, "   <-- outside GAUSS bound [-100,100]", "")
display as text "  sigma  = " %12.5f `gsig' cond(`gsig'/${JD_sc} > 100, "   <-- outside GAUSS bound [0,100]", cond(`gsig'/${JD_sc} < 0.05, "   <-- collapsing toward 0: likely degenerate likelihood", ""))
display as text "  lambda = " %12.5f `glam' cond(`glam' < 1e-3, "   <-- at 0: theta and del are then not identified", cond(`glam' > 1-1e-3, "   <-- at upper bound 1", ""))
display as text "  theta  = " %12.5f `gth'  cond(abs(`gth'/${JD_sc}) >= 100, "   <-- outside GAUSS bound [-100,100]", "")
display as text "  del    = " %12.5f `gdel' cond(`gdel'/${JD_sc} > 100, "   <-- outside GAUSS bound [0,100]", "")

* natural-scale estimates (delta-method SEs), ready for a table
estimates restore raw_jd
capture noisily nlcom (Mu: _b[mu:_cons]) (Sigma: exp(_b[/lnsigma]))       ///
    (Lambda: invlogit(_b[/lamlgt])) (Theta: _b[/theta]) (Del: exp(_b[/lndel])), post
ereturn scalar ll_ours = ll_jd
ereturn scalar N_ll    = n_jd
ereturn scalar k_par   = k_jd
* >>>>>>>>>> ESTIMATES STORE (placeholder 2 of 7): jump-diffusion, no GARCH <<<<<<<<<<
estimates store jd

}

*===============================================================================
* SECTION 2.  PLAIN GAUSSIAN GARCH(1,1), no jumps.  GAUSS: solar_garch (proc garch).
*===============================================================================
if "${SC_run_gg}" == "1" {
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

estimates store raw_garch
scalar ll_garch = e(ll)
scalar k_garch  = cond("${GG_fix}" == "", 4, 3)
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
scalar ll_arch  = .
scalar ll_archt = .
scalar df_archt = .
scalar a_archt  = .
scalar b_archt  = .
scalar lyap_t   = .
if "${GG_benchmark}" == "1" {
    display as text _n "Benchmark: Stata arch GARCH(1,1) (all rows, NO persistence cap)"
    capture noisily arch oilp, arch(1) garch(1) nolog
    if _rc == 0 {
        scalar ll_arch = e(ll)
        scalar n_arch  = e(N)
        scalar k_arch  = 4
        * >>>>>>>>>> ESTIMATES STORE (placeholder 3 of 7): Stata arch, Gaussian GARCH(1,1) <<<<<<<<<<
        estimates store arch_gauss
        display as text "*** The log likelihood in the table above is ARCH's, not ours.  Ours is on the SUMMARY line at the end. ***"
        matrix ba = e(b)
        display as text "arch coefficient vector (const, ARCH L1, GARCH L1, variance const):"
        matrix list ba, format(%12.6f) noheader
        display as text "(ours: mu, then below kappa/beta/alpha; arch's ARCH term = our alpha, its GARCH term = our beta, its variance constant = our kappa)"
    }
}

*-------------------------------------------------------------------------------
* Benchmark: GARCH(1,1) with Student-t innovations (arch, distribution(t)).
* Heavy tails WITHOUT jumps.  If persistence falls below 1 here too, fat-tailed
* errors alone remove the explosive behaviour and the jump model must be
* justified against this alternative, not only against the Gaussian GARCH.
*-------------------------------------------------------------------------------
if "${GG_benchmark_t}" == "1" {
    display as text _n "Benchmark: Stata arch GARCH(1,1) with Student-t innovations (all rows, NO persistence cap)"
    capture noisily arch oilp, arch(1) garch(1) distribution(t) nolog
    if _rc == 0 {
        scalar ll_archt = e(ll)
        scalar n_archt  = e(N)
        scalar k_archt  = 5
        * >>>>>>>>>> ESTIMATES STORE (placeholder 4 of 7): Stata arch, Student-t GARCH(1,1) <<<<<<<<<<
        estimates store arch_t
        matrix bt = e(b)
        local cnt : colfullnames bt
        display as text "arch(t) coefficient names: `cnt'"
        scalar a_archt = bt[1,2]               // ARCH term  (assumed layout: mean const, ARCH, GARCH, variance const, shape)
        scalar b_archt = bt[1,3]               // GARCH term
        capture scalar df_archt = e(tdf)
        if df_archt == . {
            display as text "(could not read the t degrees of freedom from e(tdf): read it from the table above)"
        }
        mata: st_numscalar("lyap_t", gg_lyap_t(st_numscalar("b_archt"), st_numscalar("a_archt"), st_numscalar("df_archt")))
        display as text _n "GARCH-t: ll = " %12.5f ll_archt "   df = " %9.4f df_archt
        display as text "         alpha = " %9.6f a_archt "   beta = " %9.6f b_archt "   alpha+beta = " %9.6f (a_archt + b_archt)
        display as text "         Nelson exponent with t innovations = " %9.5f lyap_t "   (> 0: not strictly stationary)"
        if (a_archt + b_archt) < 1 & ll_archt < . {
            display as text "         Persistence is below 1 with t errors: heavy tails alone remove the explosive behaviour."
        }
        else if ll_archt < . {
            display as text "         Persistence is still >= 1 with t errors: heavy tails alone do not remove it."
        }
        display as text "(Check the coefficient names above: the positions of ARCH and GARCH terms are assumed from arch's usual layout.)"
    }
}
estimates restore raw_garch

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
      (Persist: `S'), post
ereturn scalar ll_ours = ll_garch
ereturn scalar N_ll    = n_garch
ereturn scalar k_par   = k_garch
* >>>>>>>>>> ESTIMATES STORE (placeholder 5 of 7): Gaussian GARCH(1,1), our GAUSS-faithful likelihood <<<<<<<<<<
estimates store garch

mata: gg_bound_report(st_matrix("braw"))

display as text _n "Units: data scaled by " ${GG_scale} ".  To compare with estimates in other units: mu scales by the ratio,"
display as text "kappa and the unconditional variance by its square; beta and alpha are unit-free."

*-------------------------------------------------------------------------------
* One unambiguous line per run: copy this when comparing runs
*-------------------------------------------------------------------------------
mata: st_matrix("bsum", gg_natural(st_matrix("braw")))
mata: st_numscalar("gg_lyap", gg_lyap(st_matrix("bsum")[1,3], st_matrix("bsum")[1,4]))
display as text _n "SUMMARY-T | arch with Student-t: ll=" %12.5f ll_archt " df=" %9.4f df_archt " alpha=" %9.6f a_archt        ///
    " beta=" %9.6f b_archt " sum=" %9.6f (a_archt + b_archt) " nelson=" %9.5f lyap_t
display as text _n "SUMMARY | cap=" %6.4f gg_cap " | GG_fix=${GG_fix} | rc=" gg_rc                          ///
    " | OURS: ll=" %12.5f ll_garch " beta=" %9.6f bsum[1,3] " alpha=" %9.6f bsum[1,4]                    ///
    " sum=" %9.6f (bsum[1,3] + bsum[1,4]) " kappa=" %12.8f bsum[1,2] " nelson=" %8.5f gg_lyap " | ARCH benchmark ll=" %12.5f ll_arch

* Append the same facts to a file, so a sequence of runs can be read off one place
* (written next to your do-file's working directory; delete the file to start fresh).
capture {
    tempname fh
    file open `fh' using "solar_runs.txt", write append text
    file write `fh' "`c(current_date)' `c(current_time)' | cap=" %6.4f (gg_cap) " | GG_fix=${GG_fix} | rc=" (gg_rc)  ///
        " | ours ll=" %13.5f (ll_garch) " beta=" %10.7f (bsum[1,3]) " alpha=" %10.7f (bsum[1,4])                    ///
        " sum=" %10.7f (bsum[1,3] + bsum[1,4]) " kappa=" %13.9f (bsum[1,2]) " nelson=" %9.5f (gg_lyap) " | arch ll=" %13.5f (ll_arch)           ///
        " | archT ll=" %13.5f (ll_archt) " df=" %9.4f (df_archt) " alpha=" %9.6f (a_archt) " beta=" %9.6f (b_archt) " nelsonT=" %9.5f (lyap_t) ///
        " | scale=${GG_scale} droptails=${SC_droptails} dropzero=${SC_dropzero} keepfirst=${SC_keepfirst}" _n
    file close `fh'
}
display as text "Run appended to solar_runs.txt in " c(pwd)

}

*-------------------------------------------------------------------------------
* GARCH once more with the persistence cap effectively removed.  The capped fit
* (GG_cap, default .99) is the coauthor's GAUSS baseline; the uncapped fit is the
* best Gaussian GARCH and the fair comparator for the jump model.
*-------------------------------------------------------------------------------
if "${SC_run_gg}" == "1" & "${GG_extra_uncapped}" == "1" {
    display as text _n "GARCH(1,1) again with the cap removed (GG_cap_u = ${GG_cap_u}, GG_fix empty)"
    global GG_fix_save "${GG_fix}"
    global GG_cap_save "${GG_cap}"
    global GG_fix ""
    global GG_cap ${GG_cap_u}
    scalar gg_cap = ${GG_cap}

    gg_startvec ${GG_mu0} ${GG_kappa0} ${GG_beta0} ${GG_alpha0}
    gg_fit "${GG_maxopts}"
    if gg_rc != 0 {
        display as error "WARNING: the uncapped GARCH fit did not converge (rc = " gg_rc ")."
    }
    estimates store raw_garch_u
    scalar ll_garch_u = e(ll)
    scalar n_garch_u  = e(N) - 1
    scalar k_garch_u  = 4
    matrix braw = e(b)

    local km  = gg_kmin
    local cap = gg_cap
    local S "(`cap'*invlogit(_b[/apers]))"
    local W "invlogit(_b[/ashare])"
    capture noisily nlcom (Mu:      _b[mu:_cons])                           ///
          (Kappa:   `km' + exp(_b[/lnkappa]))                               ///
          (Beta:    `S'*`W')                                                ///
          (Alpha:   `S'*(1 - `W'))                                          ///
          (Persist: `S'), post
    ereturn scalar ll_ours = ll_garch_u
    ereturn scalar N_ll    = n_garch_u
    ereturn scalar k_par   = k_garch_u
    * >>>>>>>>>> ESTIMATES STORE (placeholder 6 of 7): Gaussian GARCH(1,1), cap removed <<<<<<<<<<
    estimates store garch_u

    mata: gg_bound_report(st_matrix("braw"))
    display as text "Uncapped GARCH log likelihood: " %12.5f ll_garch_u "   (capped: " %12.5f ll_garch ")"

    global GG_fix "${GG_fix_save}"
    global GG_cap "${GG_cap_save}"
    scalar gg_cap = ${GG_cap}
}

*===============================================================================
* SECTION 3.  GARCH(1,1) WITH POISSON JUMPS.  GAUSS: solar_garch_jump (proc garchmj).
*===============================================================================
if "${SC_run_gj}" == "1" {
display as text _n "GARCH(1,1) with Jumps estimation"

* GAUSS start vector: Mu | Kappa | Beta | Alpha | Lambda | Theta | Del
gj_startvec ${GJ_mu0} ${GJ_kappa0} ${GJ_beta0} ${GJ_alpha0} ${GJ_lam0} ${GJ_theta0} ${GJ_del0}

if "${GJ_multistart}" == "1" {
    * Mixture likelihoods are multimodal: screen a grid of starts, keep the best.
    scalar gj_best = -1e300
    matrix bn = J(1, 7, .)
    foreach th in -30 -5 0.5 5 30 {
        foreach dl in 3.5 60 {
            foreach bt in 0.05 0.6 {
                local al = cond(`bt' > 0.3, 0.2, 0.3)
                gj_startvec ${GJ_mu0} ${GJ_kappa0} `bt' `al' ${GJ_lam0} `th' `dl'
                gj_fit "difficult iterate(60)" quiet
                local ll = e(ll)
                capture mata: st_matrix("bn", gj_natural(st_matrix("e(b)")))
                display as text "start theta=`th' del=`dl' beta=`bt':  ll = " %11.4f `ll'          ///
                    "  kappa=" %9.5f bn[1,2] "  beta=" %6.4f bn[1,3] "  alpha=" %6.4f bn[1,4]    ///
                    "  lambda=" %6.4f bn[1,5] "  theta=" %8.4f bn[1,6] "  del=" %8.4f bn[1,7]   ///
                    "  rc=" gj_rc
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

gj_fit "${GJ_maxopts}"
if gj_rc != 0 {
    display as error _n "WARNING: ml did not converge (rc = " gj_rc ").  Treat the results below with suspicion;"
    display as error "the bound report at the end shows whether a constraint is active."
}

estimates store raw_garchjump
scalar ll_jump = e(ll)
scalar k_jump  = 7
scalar n_jump  = e(N) - 1                  // the seed row contributes nothing
matrix braw = e(b)                       // raw (unconstrained) estimates, kept for the bound check

display as text _n "N in the likelihood = " n_jump " (ml reports " (n_jump+1) " rows incl. the seed row);  log likelihood = " %11.4f ll_jump "  (data units as scaled)"
display as text "Log likelihood in raw-return units (comparable across GJ_scale): " %11.4f (ll_jump + n_jump*ln(${GJ_scale}))

if ll_garch < . {
    display as text _n "Gaussian GARCH (section 2, our GAUSS-faithful likelihood): " %11.4f ll_garch "   (cap = ${GG_cap}, GG_fix = ${GG_fix})"
    display as text   "Jump-GARCH minus GARCH:                  " %11.4f (ll_jump - ll_garch) "  (3 extra parameters: lambda, theta, del)"
    if ll_jump < ll_garch - 1 {
        display as error "Jump model is WORSE than its own special case: a local optimum or a bug.  Do not use these estimates."
    }
}
if ll_garch_u < . {
    display as text   "Gaussian GARCH, cap removed:             " %11.4f ll_garch_u
    display as text   "Jump-GARCH minus uncapped GARCH:         " %11.4f (ll_jump - ll_garch_u) "  <-- the fair comparison"
}
if ll_garch < . | ll_garch_u < . {
    display as text "(Not a chi2 test: at lambda = 0, theta and del are unidentified.)"
}

*-------------------------------------------------------------------------------
* Self-check: recompute the log likelihood at the estimates with a plain Stata
* loop (no Mata, no log-sum-exp) and compare with ml's value.
*-------------------------------------------------------------------------------
if "${GJ_selfcheck}" == "1" {
    mata: st_matrix("bnat", gj_natural(st_matrix("braw")))
    local cmu  = bnat[1,1]
    local ckap = bnat[1,2]
    local cbet = bnat[1,3]
    local calp = bnat[1,4]
    local clam = bnat[1,5]
    local cth  = bnat[1,6]
    local cdel = bnat[1,7]

    * (a) plain-Stata recursion and direct (non-log) mixture formula
    quietly {
        generate double _chk_u2 = (oilp - `cmu')^2
        summarize _chk_u2
        generate double _chk_h = r(mean) in 1
        forvalues i = 2/`=_N' {
            replace _chk_h = `ckap' + `calp'*_chk_u2[`i'-1] + `cbet'*_chk_h[`i'-1] in `i'
        }
        generate double _chk_d = 0
        forvalues k = 0/$GJ_K {
            local fk = exp(lnfactorial(`k'))        // k!  (Stata has no factorial() function)
            replace _chk_d = _chk_d + (`clam'^`k'/`fk') * (_chk_h + `k'*`cdel'^2)^(-0.5) ///
                * exp(-0.5*(oilp - `cmu' - `k'*`cth')^2/(_chk_h + `k'*`cdel'^2))
        }
        generate double _chk_ll = -`clam' - 0.5*ln(2*_pi) + ln(_chk_d)
        summarize _chk_ll if t > 1
    }
    local llchk = r(sum)
    local nchk  = r(N)

    * (b) the Mata likelihood called directly at the same raw estimates
    tempname q1 q2 q3 q4 q5 q6 q7
    forvalues i = 1/7 {
        scalar `q`i'' = braw[1,`i']
    }
    quietly {
        generate double _chk_lm = .
        mata: gj_ll("_chk_lm", "oilp", ("`q1'","`q2'","`q3'","`q4'","`q5'","`q6'","`q7'"), $GJ_K)
        summarize _chk_lm if t > 1
    }
    local llm  = r(sum)
    local nllm = r(N)

    display as text _n "Self-check, three numbers that should agree:"
    display as text "  ml e(ll)                          = " %14.6f ll_jump "   (N = " n_jump ")"
    display as text "  Mata gj_ll summed over t>1        = " %14.6f `llm'  "   (N = `nllm')"
    display as text "  plain-Stata loop summed over t>1  = " %14.6f `llchk' "   (N = `nchk')"

    quietly summarize _chk_lm
    local llall = r(sum)
    display as text "  Mata gj_ll summed over ALL rows   = " %14.6f `llall' "   (N = " r(N) ")"
    display as text "  Mata lnf at t=1: " %12.6f _chk_lm[1] "   at t=2: " %12.6f _chk_lm[2] "   at t=N: " %12.6f _chk_lm[_N]
    display as text "  e(ll) minus Mata sum over t>1     = " %12.6f (ll_jump - `llm')

    quietly {
        generate double _chk_dif = _chk_lm - _chk_ll
        count if missing(_chk_ll) & t > 1
        local nmiss = r(N)
        summarize _chk_dif if t > 1
        local dmean = r(mean)
        local dmin  = r(min)
        local dmax  = r(max)
        summarize t if abs(_chk_dif) > 1e-8 & t > 1
        local ndif = r(N)
        local tfirst = r(min)
        local tlast  = r(max)
    }
    display as text "  rows where the plain loop is missing: `nmiss'"
    display as text "  per-observation Mata minus loop: mean " %12.8f `dmean' "  min " %12.8f `dmin' "  max " %12.8f `dmax'
    display as text "  rows differing by more than 1e-8: `ndif'" cond(`ndif' > 0, " (first at t = `tfirst', last at t = `tlast')", "")
    if `ndif' > 0 {
        display as text "  the ten largest differences:"
        gsort -_chk_dif
        list t oilp _chk_lm _chk_ll _chk_dif in 1/5, noobs
        gsort _chk_dif
        list t oilp _chk_lm _chk_ll _chk_dif in 1/5, noobs
        sort t
    }

    if abs(ll_jump - `llm') > 1e-6*max(1, abs(ll_jump)) {
        display as error "ml's e(ll) differs from the sum of the Mata likelihood: an ml sample/weighting issue, not a formula issue."
    }
    if abs(`llm' - `llchk') > 1e-6*max(1, abs(`llm')) {
        display as error "SELF-CHECK FAILED: Mata likelihood and plain-Stata recomputation disagree (see the row pattern above)."
    }
    else display as text "Self-check passed."
    drop _chk_*
}

*-------------------------------------------------------------------------------
* Results on the original (GAUSS) parameter scale, delta-method SEs
* (GAUSS: print b; print sqrt(diag(h)))
*-------------------------------------------------------------------------------
local S "(2*0.0001 + (0.99 - 2*0.0001)*invlogit(_b[/apers]))"
local W "invlogit(_b[/ashare])"

estimates restore raw_garchjump
local cn : colfullnames e(b)
display as text _n "coefficient names: `cn'"

capture noisily nlcom (Mu:      _b[mu:_cons])                                            ///
      (Kappa:   exp(_b[/lnkappa]))                                  ///
      (Beta:    0.0001 + (`S' - 2*0.0001)*`W')                           ///
      (Alpha:   0.0001 + (`S' - 2*0.0001)*(1 - `W'))                     ///
      (Lambda:  0.0001 + (1 - 0.0001)*invlogit(_b[/lamlgt]))        ///
      (Theta:   _b[/theta])                                         ///
      (Del:     0.0001 + exp(_b[/lndel]))                           ///
      (Persist: `S'), post
ereturn scalar ll_ours = ll_jump
ereturn scalar N_ll    = n_jump
ereturn scalar k_par   = k_jump
* >>>>>>>>>> ESTIMATES STORE (placeholder 7 of 7): GARCH(1,1) with Poisson jumps <<<<<<<<<<
estimates store garchjump

mata: gj_bound_report(st_matrix("braw"))

display as text _n "Units: data scaled by " ${GJ_scale} ".  To compare with estimates in other units: mu, theta, del scale by the ratio,"
display as text "kappa by its square; beta, alpha, lambda are unit-free."

}

*===============================================================================
* MODEL COMPARISON
* Log likelihoods in raw-return units (ll + N*ln(scale)), so runs at different
* SC_scale are comparable.  N = observations actually in the likelihood.
* The GARCH-family rows (garch, garch_u, jump) use N-1 observations (seed row);
* gbm/jd use N unless JD_skipseed = 1; the arch rows use all N rows and their own
* variance initialisation.  Compare log likelihoods across rows only when the
* samples agree, and read AIC/BIC as rough guides, not tests (Poisson mixtures
* violate the usual regularity conditions).
*===============================================================================
display as text _n "{hline 98}"
display as text "Model comparison (log likelihood in raw-return units)"
display as text "{hline 98}"
display as text %-30s "model" %6s "k" %8s "N" %14s "loglik" %14s "AIC" %14s "BIC"
foreach m in gbm jd garch garch_u jump arch archt {
    local lab "`m'"
    if "`m'" == "gbm"     local lab "GBM (normal)"
    if "`m'" == "jd"      local lab "Jump-diffusion"
    if "`m'" == "garch"   local lab "GARCH, cap = ${GG_cap}"
    if "`m'" == "garch_u" local lab "GARCH, cap removed"
    if "`m'" == "jump"    local lab "GARCH + jumps"
    if "`m'" == "arch"    local lab "arch Gaussian (Stata)"
    if "`m'" == "archt"   local lab "arch Student-t (Stata)"
    local lr  = scalar(ll_`m') + scalar(n_`m')*ln(${SC_scale})
    local aic = -2*`lr' + 2*scalar(k_`m')
    local bic = -2*`lr' + scalar(k_`m')*ln(scalar(n_`m'))
    display as text %-30s "`lab'" %6.0f scalar(k_`m') %8.0f scalar(n_`m') %14.4f `lr' %14.4f `aic' %14.4f `bic'
}
display as text "{hline 98}"

*===============================================================================
* TABLE TEMPLATES.  The stored estimates are (natural scale, delta-method SEs):
*     gbm   jd   garch   garch_u   garchjump      and, from Stata's arch,  arch_gauss   arch_t
* Parameter names line up across models (Mu, Sigma, Kappa, Beta, Alpha, Lambda, Theta, Del, Persist).
* Uncomment ONE of these.  They have not been run; adjust options to taste.
*===============================================================================
* estimates table gbm jd garch garch_u garchjump, b(%10.5f) se stats(ll_ours N_ll k_par)
*
* etable, estimates(gbm jd garch garch_u garchjump) mstat(N_ll) mstat(ll_ours) mstat(k_par) ///
*     showstars showstarsnote title("Solar SREC models") export(solar_models.docx, replace)
*
* ssc install estout          // once
* esttab gbm jd garch garch_u garchjump using solar_models.rtf, replace se b(%9.4f) ///
*     scalars("ll_ours Log likelihood" "N_ll Obs. in likelihood" "k_par Parameters") ///
*     mtitles("GBM" "Jump-diffusion" "GARCH" "GARCH (no cap)" "GARCH+jumps")

display as text _n "solar_all.do finished."
