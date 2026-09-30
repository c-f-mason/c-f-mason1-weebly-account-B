/*==============================================================================
  gj_probe.do
  Diagnostic.  Run in the SAME Stata session, right after solar_garch_jump.do
  has finished (it needs the data in memory, the matrix braw, and the Mata
  functions gj_natural / gj_ll).

  Why: ml evaluates the log likelihood at braw as 200.08043, but calling the
  same Mata function directly at braw gives 199.82041 (row values: lnf[1] =
  -0.25490772, lnf[2] = 0.26622989, lnf[2101] = 1.2512374).  This runs the
  evaluator under ml and prints what it sees, so we can see where the two differ.

  Run with:   do gj_probe.do
==============================================================================*/

version 19

display as text _n "braw (the parameters being tested), full precision:"
matrix list braw, format(%18.12f)

capture program drop gj_probe
program define gj_probe
    version 19
    args todo b lnfj
    tempname p1 p2 p3 p4 p5 p6 p7
    forvalues i = 1/7 {
        mleval `p`i'' = `b', eq(`i') scalar
    }
    display as text _n "probe: _N = " _N "   t[1] = " t[1] "   t[2] = " t[2] "   t[_N] = " t[_N]
    display as text "probe: oilp[1] = " %14.10f oilp[1] "   oilp[2] = " %14.10f oilp[2] "   oilp[_N] = " %14.10f oilp[_N]
    display as text "probe: params received from ml:"
    display as text "   " %16.10f `p1' %16.10f `p2' %16.10f `p3' %16.10f `p4' %16.10f `p5' %16.10f `p6' %16.10f `p7'
    mata: gj_ll("`lnfj'", "oilp", ("`p1'","`p2'","`p3'","`p4'","`p5'","`p6'","`p7'"), 10)
    quietly summarize `lnfj' if $ML_samp
    display as text "probe: sum over ML sample  = " %16.8f r(sum) "   N = " r(N)
    quietly summarize `lnfj'
    display as text "probe: sum over all rows   = " %16.8f r(sum) "   N = " r(N)
    display as text "probe: lnfj[1] = " %14.8f `lnfj'[1] "   lnfj[2] = " %14.8f `lnfj'[2] "   lnfj[_N] = " %14.8f `lnfj'[_N]
    display as text "probe: (direct call gave lnf[1] = -0.25490772, lnf[2] = 0.26622989, lnf[2101] = 1.2512374)"
end

ml model gf0 gj_probe (mu: oilp = ) /lnkappa /apers /ashare /lamlgt /theta /lndel if t > 1
ml init braw, copy
capture noisily ml maximize, iterate(0)
