/*==============================================================================
  solar_garch_jump.do
  Stata 19 translation of the GAUSS program "solar_garch_jump" (Solar data
  analysis, 02/22/17; data switch updated May 2023).

  Model (the likelihood actually optimised in GAUSS is proc garchmj):

      y_t | K_t = k  ~  N( mu + k*theta ,  h_t + k*del^2 ),   K_t ~ Poisson(lambda)
      h_t = kappa + alpha*(y_{t-1}-mu)^2 + beta*h_{t-1},      h_1 = sample mean of (y-mu)^2

  i.e. a GARCH(1,1) diffusion with a constant-intensity Poisson jump component,
  the Poisson sum truncated at K = 10 jumps per day.  Parameter order and names
  match the GAUSS program: Mu Kappa Beta Alpha Lambda Theta Del.

  Not carried over (dead code in the GAUSS file): proc garch (single-jump
  Bernoulli version), proc grd (analytic GARCH gradient, its use is commented
  out), _dshape, and the unused h0..h10 recursions inside garchmj.

  Lessons carried over from validating the no-GARCH jump program (solar_jump.do)
  against the coauthor's output:
   * SCALE.  Estimates are scale-equivariant: if the data are multiplied by c,
     mu, theta, del scale by c, kappa by c^2, and beta, alpha, lambda do not
     move.  GJ_scale sets the data units; start values and the screening grid
     are written in GJ_startscale units and rescaled automatically.  Compare
     with the coauthor only after converting units.  The log likelihood shifts
     by N*ln(c); the script prints it in raw-return units.
   * SAMPLE.  Sample mismatches (2101 vs 2099 rows) were the main obstacle.
     GJ_dropzero, GJ_droptails, GJ_keepfirst reproduce the ways a sample can
     differ (see settings).  Note for GARCH: deleting rows from the middle of
     the series (dropzero) makes non-adjacent days adjacent in the recursion.
   * ZEROS.  Exact zeros are listed.  Unlike the constant-variance model,
     kappa>0 keeps h_t away from 0 here, so the sigma->0 spike is much weaker,
     but watch kappa in the screening table.
   * MULTIMODALITY.  Single-start runs stalled at poor optima (-9583, -9589
     vs -9470.78 after screening).  Screening is on by default.
   * BENCHMARK.  Plain GARCH(1,1) via Stata's arch is fitted first; the jump
     model nests it (lambda -> 0), so its log likelihood must not be lower.
   * SELF-CHECK.  There is no GAUSS GARCH output to compare against, so the
     Mata likelihood is re-computed at the estimates with a plain Stata loop
     and the two are compared.

  Differences from GAUSS/CML you should know about
   * CML enforces bounds and beta+alpha<=.99 directly.  Stata's ml cannot, so
     the constraints are imposed by smooth reparameterisation (see gj_natural()
     below).  Results are identical when the optimum is interior; if a
     constraint binds, the transformed parameter runs off to +/-infinity and
     the delta-method SEs are unreliable.  A bound check is printed at the end.
   * CML's Lagrange multiplier printout has no counterpart here.
   * GAUSS zero-weights obs 1 and rescales the others by N/(N-1); here obs 1 is
     simply excluded (if t>1).  Point estimates are unaffected; SEs differ by a
     factor of about sqrt((N-1)/N).
   * The Poisson mixture is summed in logs (log-sum-exp) for numerical
     stability.  Mathematically the same as the GAUSS expression.
   * Not yet compared with GAUSS GARCH output.  Obtain the coauthor's Model
     estimates, convert units, and match the sample (Model 1 of solar_jump.do
     fingerprints a sample) before relying on it.

  Open modelling question (unchanged): h_t is driven by the RAW squared residual
  (y-mu)^2, so a jump day feeds into next-period variance at full size.
==============================================================================*/

version 19
clear all
macro drop GJ_*          // start from clean settings; stale globals from an earlier run would otherwise persist
set more off

*-------------------------------------------------------------------------------
* User settings
*-------------------------------------------------------------------------------
global GJ_datafile "/Users/chuckmason/Dropbox/Research/NeilWilmot_jumps/SREC/PJM_SREC_prices.dta"
global GJ_var      SRECp_ret      // variable to analyse
global GJ_scale    1              // oilp = GJ_scale * GJ_var.  1 = raw (matches the jump-only coauthor output), 100 = percent
global GJ_startscale 100          // data units in which the start values / screening grid below are written (GAUSS: 10 for
                                  //   solar=4 residuals; the screening grid was tuned at 100)
global GJ_sortvar  ""             // optional: date variable to sort by; the GARCH recursion needs time order

* sample switches, applied in this order
global GJ_droptails  0            // 1 = drop first and last USABLE observations (after missing returns are removed)
                                  // 2 = drop first and last ROWS of the file (before missing returns are removed)
global GJ_dropzero   0            // 1 = drop rows with an exact zero return
global GJ_keepfirst  0            // n>0 = keep only the first n rows, as GAUSS's load solmat[n,k] does

// options passed to ml maximize.  Defaults: tolerance(1e-6) ltolerance(1e-7) nrtolerance(1e-5)
global GJ_maxopts "difficult iterate(100) showtolerance"
global GJ_multistart 1         // 1 = screen a grid of start values first, then polish the best (GAUSS used one start)
global GJ_benchmark  1         // 1 = also fit plain GARCH(1,1) with -arch- and compare
global GJ_selfcheck  1         // 1 = recompute the log likelihood with a plain Stata loop and compare
global GJ_K   10               // maximum number of jumps per period in the Poisson sum

* start values, GAUSS order: Mu | Kappa | Beta | Alpha | Lambda | Theta | Del   (in GJ_startscale units)
global GJ_mu0 = 0.10
global GJ_kappa0 = 1.5
global GJ_beta0 = 0.05
global GJ_alpha0 = 0.3
global GJ_lam0 = 0.3
global GJ_theta0 = 0.50
global GJ_del0 = 3.5

*-------------------------------------------------------------------------------
* Data (loaded directly)
*-------------------------------------------------------------------------------
display as text "Solar Jump paper estimation"
display as text " GARCH(1,1) with Poisson jumps"

use "${GJ_datafile}", clear
if "${GJ_sortvar}" != "" sort ${GJ_sortvar}
gen double oilp = ${GJ_scale}*${GJ_var}
if "${GJ_droptails}" == "2" {
    drop in 1
    drop in l
}
drop if oilp == .              // the recursion needs an unbroken series
if "${GJ_droptails}" == "1" {
    drop in 1
    drop in l
}
if "${GJ_dropzero}" == "1" drop if oilp == 0
if ${GJ_keepfirst} > 0 keep in 1/${GJ_keepfirst}
display as text "Observations used: " _N

gen long t = _n                // row order is time order (check this!)
tsset t

summarize oilp, detail

quietly count if oilp == 0
local nzero = r(N)
display as text _n "Exact zero returns: `nzero' of " _N
if `nzero' > 0 {
    display as text "Rows with an exact zero return (unchanged price, or a missing value coded as 0?):"
    list t ${GJ_var} if oilp == 0, noobs
}

global GJ_Y oilp
global GJ_sc = ${GJ_scale}/${GJ_startscale}     // start values are multiplied by this (kappa by its square)

*-------------------------------------------------------------------------------
* Mata: parameter transformation and log-likelihood
*-------------------------------------------------------------------------------
* Unconstrained -> constrained (GAUSS bounds):
*   kappa  > 0                          kappa = exp(r2)
*   beta, alpha >= 1e-4, beta+alpha <= .99   (s = beta+alpha, w = beta's share)
*   1e-4 <= lambda <= 1
*   del >= 1e-4
*   mu, theta unrestricted
mata:
mata clear

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

*-------------------------------------------------------------------------------
* ml evaluator (method gf0: observation-level log likelihood, recursive in t)
*-------------------------------------------------------------------------------
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

*-------------------------------------------------------------------------------
* Guard against a partial run: if the settings block at the top was skipped,
* stop with a message instead of expanding empty macros.
*-------------------------------------------------------------------------------
if "${GJ_K}" == ""      global GJ_K 10
if "${GJ_Y}" == ""      global GJ_Y oilp
if "${GJ_mu0}" == ""    display as error "GJ_* start values are not set: run the WHOLE do-file (do solar_garch_jump.do), not a selection"
if "${GJ_mu0}" == ""    exit 198

*-------------------------------------------------------------------------------
* Helpers: start vector on the unconstrained scale, and one ml fit
*-------------------------------------------------------------------------------
* Arguments are in GJ_startscale units; they are converted to the data's units here.
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

* One ml fit from matrix b0.  ml maximize aborts the do-file (r(430)) when it
* does not converge; capture so the script can carry on.  Sets scalar gj_rc.
capture program drop gj_fit
program define gj_fit
    version 19
    args maxopts quiet
    ml model gf0 gj_eval (mu: oilp = ) /lnkappa /apers /ashare /lamlgt /theta /lndel ///
        if t > 1, title("GARCH(1,1) with Poisson jumps")
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

*-------------------------------------------------------------------------------
* Benchmark: plain Gaussian GARCH(1,1) with Stata's arch.  The jump model nests
* it (lambda -> 0), so the jump log likelihood should not be lower.  arch
* initialises the variance differently from the GAUSS recursion, so expect
* small differences, not identity.
*-------------------------------------------------------------------------------
scalar ll_arch = .
if "${GJ_benchmark}" == "1" {
    display as text _n "Benchmark: Gaussian GARCH(1,1) via arch (no jumps)"
    capture noisily arch oilp if t > 1, arch(1) garch(1) nolog
    if _rc == 0 {
        scalar ll_arch = e(ll)
        estimates store benchmark
    }
}

*-------------------------------------------------------------------------------
* Estimation.  GAUSS: CML, BHHH algorithm, step halving.  Stata's gf0 evaluator
* rejects technique(), so ml's default Newton-Raphson (numerical derivatives) is
* used, with 'difficult' to help in flat regions.  Same optimum, different path.
* First observation is excluded (GAUSS weight = 0) but still seeds the recursion.
*-------------------------------------------------------------------------------
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

estimates store garchjump
scalar ll_jump = e(ll)
scalar n_jump  = e(N)
matrix braw = e(b)                       // raw (unconstrained) estimates, kept for the bound check

display as text _n "N = " n_jump ";  log likelihood = " %11.4f ll_jump "  (data units as scaled)"
display as text "Log likelihood in raw-return units (comparable across GJ_scale): " %11.4f (ll_jump + n_jump*ln(${GJ_scale}))

if ll_arch < . {
    display as text _n "Benchmark arch GARCH(1,1) log likelihood: " %11.4f ll_arch
    display as text   "Jump-GARCH minus GARCH:                  " %11.4f (ll_jump - ll_arch) "  (3 extra parameters: lambda, theta, del)"
    if ll_jump < ll_arch - 1 {
        display as error "Jump model is WORSE than its own special case: a local optimum or a bug.  Do not use these estimates."
    }
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

    quietly {
        generate double _chk_u2 = (oilp - `cmu')^2
        summarize _chk_u2
        generate double _chk_h = r(mean) in 1
        forvalues i = 2/`=_N' {
            replace _chk_h = `ckap' + `calp'*_chk_u2[`i'-1] + `cbet'*_chk_h[`i'-1] in `i'
        }
        generate double _chk_d = 0
        forvalues k = 0/$GJ_K {
            replace _chk_d = _chk_d + (`clam'^`k'/factorial(`k')) * (_chk_h + `k'*`cdel'^2)^(-0.5) ///
                * exp(-0.5*(oilp - `cmu' - `k'*`cth')^2/(_chk_h + `k'*`cdel'^2))
        }
        generate double _chk_ll = -`clam' - 0.5*ln(2*_pi) + ln(_chk_d)
        summarize _chk_ll if t > 1
    }
    local llchk = r(sum)
    display as text _n "Self-check: Mata/ml log likelihood = " %14.6f ll_jump "   plain-Stata loop = " %14.6f `llchk'
    if abs(ll_jump - `llchk') > 1e-6*max(1, abs(ll_jump)) {
        display as error "SELF-CHECK FAILED: the Mata likelihood and the plain-Stata recomputation disagree."
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

mata: gj_bound_report(st_matrix("braw"))

display as text _n "Units: data scaled by " ${GJ_scale} ".  To compare with estimates in other units: mu, theta, del scale by the ratio,"
display as text "kappa by its square; beta, alpha, lambda are unit-free."
