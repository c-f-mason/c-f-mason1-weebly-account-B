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
   * Untested against GAUSS output: compare log likelihood and estimates
     before relying on it.
==============================================================================*/

version 19
clear all
set more off

*-------------------------------------------------------------------------------
* User settings
*-------------------------------------------------------------------------------
local solar   = 4          // 1 = SREC price returns, 2 = Henry Hub, 3 = PJM Wh,
                           // 4 = residuals from Chuck (May 2023)
local datadir "."          // folder holding the .txt files (GAUSS: C:\gauss22\Neil\solar\)
global GJ_K   10           // maximum number of jumps per period in the Poisson sum

* start values, GAUSS order: Mu | Kappa | Beta | Alpha | Lambda | Theta | Del
local mu0 = 0.10
local kappa0 = 1.5
local beta0 = 0.05
local alpha0 = 0.3
local lam0 = 0.3
local theta0 = 0.50
local del0 = 3.5

*-------------------------------------------------------------------------------
* Data
*-------------------------------------------------------------------------------
display as text "Solar Jump paper estimation"
display as text " Single Jump processes"

if `solar' == 1 {
    display as text "Solar Price"
    display as text "Start data for sample Early 08/01/2009 up to 11/30/2015"
    display as text "DAILY Data"
    infile cnt price using "`datadir'/solar.txt", clear
    local nexp = 2099
    gen double oilp = 1*price
}
else if `solar' == 2 {
    display as text "Henry Hub Price returns"
    display as text "Start data for sample Early 07/31/2009 up to 11/30/2015"
    display as text "DAILY Data"
    infile cnt price rtn using "`datadir'/HH_prices.txt", clear
    local nexp = 1599
    gen double oilp = 100*rtn
}
else if `solar' == 3 {
    display as text "PJM Wh Electricity Prices"
    display as text "Start data for sample Early 07/31/2009 up to 11/30/2015"
    display as text "DAILY Data"
    infile cnt price rtn using "`datadir'/pjm_prices.txt", clear
    local nexp = 1606
    gen double oilp = 100*rtn
}
else if `solar' == 4 {
    display as text "Residuals Data from Chuck -- May 2023"
    display as text "Start data for sample Early 07/31/2009 up to 11/30/2015"
    display as text "DAILY Data"
    infile cnt resid using "`datadir'/udata.txt", clear
    local nexp = 2099
    gen double oilp = 10*resid
}
else {
    display as error "PROBLEM WITH DATA"
    exit 198
}

if _N != `nexp' {
    display as error "warning: expected `nexp' observations, found " _N
}
assert !missing(oilp)          // the recursion needs an unbroken series

gen long t = _n                // file order is time order, as in GAUSS
tsset t

summarize oilp, detail

global GJ_Y oilp

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
    printf("\n{txt}Natural-scale estimates and bound check\n")
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
        printf("{err}  A bound is active: the Hessian-based SEs above are not valid.\n")
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
* Start values on the unconstrained scale
*-------------------------------------------------------------------------------
local lb = 0.0001
local ub = 0.99
local s0 = `beta0' + `alpha0'

matrix b0 = ( `mu0',                                          ///
              ln(`kappa0'),                                   ///
              logit((`s0' - 2*`lb')/(`ub' - 2*`lb')),         ///
              logit((`beta0' - `lb')/(`s0' - 2*`lb')),        ///
              logit((`lam0' - `lb')/(1 - `lb')),              ///
              `theta0',                                       ///
              ln(`del0' - `lb') )

*-------------------------------------------------------------------------------
* Estimation.  GAUSS: CML, BHHH algorithm, step halving.
* First observation is excluded (GAUSS weight = 0) but still seeds the recursion.
*-------------------------------------------------------------------------------
display as text _n "GARCH(1,1) with Jumps estimation"

ml model gf0 gj_eval (mu: oilp = ) /lnkappa /apers /ashare /lamlgt /theta /lndel ///
    if t > 1, title("GARCH(1,1) with Poisson jumps")
ml init b0, copy
ml maximize, technique(bhhh 30 bfgs 30) vce(oim) difficult

estimates store garchjump

*-------------------------------------------------------------------------------
* Results on the original (GAUSS) parameter scale, delta-method SEs
* (GAUSS: print b; print sqrt(diag(h)))
*-------------------------------------------------------------------------------
matrix braw = e(b)                       // raw (unconstrained) estimates, kept for the bound check

local S "(2*0.0001 + (0.99 - 2*0.0001)*invlogit(_b[apers:_cons]))"
local W "invlogit(_b[ashare:_cons])"

nlcom (Mu:      _b[mu:_cons])                                            ///
      (Kappa:   exp(_b[lnkappa:_cons]))                                  ///
      (Beta:    0.0001 + (`S' - 2*0.0001)*`W')                           ///
      (Alpha:   0.0001 + (`S' - 2*0.0001)*(1 - `W'))                     ///
      (Lambda:  0.0001 + (1 - 0.0001)*invlogit(_b[lamlgt:_cons]))        ///
      (Theta:   _b[theta:_cons])                                         ///
      (Del:     0.0001 + exp(_b[lndel:_cons]))                           ///
      (Persist: `S'), post

mata: gj_bound_report(st_matrix("braw"))
