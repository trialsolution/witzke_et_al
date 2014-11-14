********************************************************************************
$ontext

   GAMS file : WITZKE_ET_AL.GMS

   @purpose  : numerical example in Witzke et al. conference paper
               "Modelling EU sugar market scenarios
               with a modified Armington approach"
   @author   :
   @date     : 12.05.14
   @since    :
   @refDoc   :
   @seeAlso  :
   @calledBy :

$offtext
********************************************************************************

$offlisting


*
*   --- The calibration approach requires two points (price/quantity) in time:
*       (1) historical observation
*       (2) expectation (assumption) on price/quantity pair
*

sets
  i           "commodity groups" /i1/
  r           "source of imports (origin)" /n, r/
  points      "calibration points" /data, assumption/
  approaches  "calibration approaches" /standard, modified/
;

alias(i, ii);
alias(r, rr);

parameters
  sigma(i)                      "substitution elasticity"
  rho(i)                        "substitution parameter"
  Y                             "income"
  mu(i, r, approaches)          "commitment parameter"
  delta(i, r, approaches)       "CES share parameter"
  U(i)                          "utility from consuming the commodity group i"
  p(i,r,points)                 "(shadow) prices"
;

table   share(i,r, points)  "budget shares"

        n.assumption            r.assumption
i1          .9                      .1
;

table   p_null(i, r)   "price in the base year"

               n          r
      i1       1         1.5
;

table   x_null(i, r)   "demand in the base year"

              n          r
      i1      100        1.E-4
*      i1      100        0
;

table    p_calib(i,r)  "hypothetical prices in the emerging trade flow (assumption)"

               n           r
      i1       1           .9
;


sigma(i) = 5;
rho(i)   =  (sigma(i) - 1) / sigma(i);
Y        = 100;

share(i, r, "data") = p_null(i,r) * x_null(i,r) / sum( (ii,rr), p_null(ii,rr) * x_null(ii,rr));


variables
  v_utility(i)                 "utility"
  v_x(i,r,points)              "import demand"
  v_delta(i,r)                 "CES share parameters (for calibration)"
  v_W(i,points)                "price index"
  v_mu(i,r)                    "commitment parameters"
  test_u(i,points)             "utility (only for testing purposes)"
;

 positive variable v_W(i,points), v_x(i,r,points);

equations
  price_index(i, points)       "price index"
  import_demand(i, r, points)  "optimal import demand, F.O.C. of expenditure minimization under fix utility"
  ces_agg(i,points)            "CES utility aggregator"
;

  price_index(i, points) $ v_W.range(i,points) ..
    v_W(i,points) =e= sum(r, v_delta(i,r) ** sigma(i) * p(i,r,points) ** (1-sigma(i)))  ** (1/(1-sigma(i)));

  import_demand(i,r,points) $ [(v_delta.range(i,r) and v_x.l(i,r,points))
                              or v_mu.range(i,r)]
                              ..
    v_x(i,r,points) =e= v_utility(i) * v_delta(i,r) ** sigma(i) * (p(i,r,points)/v_W(i,points)) ** (-sigma(i))
      + v_mu(i,r);

  ces_agg(i,points) $ test_u.range(i,points) ..

    test_u(i,points) =e= sum(r, v_delta(i,r) * (v_x(i,r,points) - v_mu(i,r))**rho(i))**(1/rho(i));


  model ces_aggregator /ces_agg/;
*
*   --- without loss of generality, one of the delta's is set to unity
*       (changing this would only cause a different scaling (monotone transformation)
*        in the utility but do not change the utility ordering)
*
*   !!! note to above
*   -----------------
*
*   it would only be the case if there is a multiplicative term (gamma)
*   in the CES aggregator. But in this version there isn't

* delta(i, "n", approaches) = 1;

*
*   ---  set benchmark utility level equal to total consumption
*        utility does not change throughout the calibration
*

 U(i)                 = sum(r, x_null(i,r));
 v_utility.fx(i)      = U(i);

*
*   ---  Calibrate the CES funciton, standard case
*        equation (5) in the paper
*
 model CES /import_demand, price_index/;

*  set observed demand
 v_x.fx(i,r,"data")          = x_null(i,r);
*  fix variables to switch of the related equations
*  (we now only calculate for the observed calibration point)
 v_x.fx(i,r,"assumption")    = 0;

*  prices
 p(i,r,"data")               = p_null(i,r);
 p(i,r,"assumption")         = 0;
 v_W.L(i,"data")             = 1;
 v_W.fx(i,"assumption")      = 1;

*  CES parameters (initialization prevents numerical problems)
 v_delta.l (i,r)            = .5;
 v_delta.lo(i,r)            = 1.E-2;
 v_delta.up(i,r)            = +inf;

*  No commitments in the standard case
 mu(i, r, "standard")        = 0;
 v_mu.fx(i,r)                = 0;

 CES.holdfixed   = 1;
 CES.solprint    = 1;
* CES.holdfixed   = 0;
* CES.iterlim     = 0;

 v_W.lo(i,points) $ v_W.range(i,points) = 1.E-2;

*
*  --- the standard calibration only works for positive demand quantities
*  --- if demand is zero then
*      we only calculate a good initial point for the modified CES calibration
*      by changing zero flows to a small poisitive value
*
if(v_x.l("i1","r","data") eq eps,

  v_x.fx("i1","r","data") = 1.E-4;
  solve CES using CNS;
  if(CES.numinfes ne 0, abort "problem with standard calibration");

  v_x.fx("i1","r","data") = 0;

else

  solve CES using CNS;
  if(CES.numinfes ne 0, abort "problem with standard calibration");

);
*$exit

* save share parameters in the standard case
if(v_x.l("i1","r","data") gt eps,

        delta(i,r,"standard") =  v_delta.L(i,r);

*  --- test if we replicate the utility level
        v_delta.fx(i,r)       =  v_delta.L(i,r);
        v_mu.fx(i,r)          =  v_mu.l(i,r);
        v_x.fx(i,r,points)    =  v_x.l(i,r,points);
        test_u.l(i,points)    =  1;
        test_u.fx(i,"assumption") = 0;
        solve ces_aggregator using CNS;


);
*abort "check standard calibration parameters", delta;
*$exit
*
*   ---  Calibrate the CES funciton, modified version
*        equation (5) in the paper
*

* -- the committment on goods from country 'n' remains zero...
  v_mu.fx(i, 'n') = 0;

* -- ...and we allow for non-zero commitment parameters for country 'r'
*    (the commitment term is negative in our case)
  v_mu.l(i, 'r')    = - 1;
  v_mu.lo(i, 'r')   = -inf;
  v_mu.up(i, 'r')   = 0;


 parameter x_calib(i,r) "demand in the second calibration point (assumption)";
 x_calib(i,r) = share(i,r,"assumption") * Y / p_calib(i,r);


  v_x.fx(i,r,"assumption")          = x_calib(i,r);
*  free the demand for region 'n' in the expected calibration point
  v_x.lo(i,"n","assumption")      = 0.1 * v_x.l(i,"n","assumption");
  v_x.up(i,"n","assumption")      = +inf;

  v_delta.lo(i,r)                 = 1.E-2;
  v_delta.up(i,r)                 = +inf;
  v_delta.l (i,r)                 = .5;

  p(i,r,"assumption")               = p_calib(i,r);

  v_W.lo(i,points)                  = 1.E-2;
  v_W.up(i,points)                  = +inf;
  v_W.l(i,points)                   = 0.5;



 CES.solprint    = 1;
 CES.holdfixed   = 0;
* CES.iterlim     = 0;


 solve CES using CNS;
 if(CES.numinfes ne 0, abort "problem with calibration model (commitment version)");
*$exit

 delta(i,r,"modified") =  v_delta.L(i,r);
 mu(i,r,"modified")    =  v_mu.L(i,r);

 option delta:4:2:1;
display "check calibration parameters", delta, mu, rho, sigma;
*$exit


*
*   ---   Test calibration [utility reproduced?]
*

        v_delta.fx(i,r)       =  v_delta.L(i,r);
        v_mu.fx(i,r)          =  v_mu.l(i,r);
        v_x.fx(i,r,points)    =  v_x.l(i,r,points);
        test_u.l(i,points)    =  1;
        test_u.lo(i,points)   =  0;
        test_u.up(i,points)   =  +inf;
        solve ces_aggregator using CNS;


*$exit

*
*   ---   Test calibration
*
*   ---   Define a separate model for simulations
*         and test it in the 'assumed' point
*
*         The reason why we need a new model is that the composite price & quantity
*         variables have an additional dimension above:
*         they are dependent on the calibration points ("points").
*
positive variables
  v_sim_W(i)       "composite price index in simulations"
  v_sim_x(i,r)     "import demand in simulations"
;

parameters
  price(i,r)  "current price in simulation (exogenous)"
;

equations
  sim_price_index(i)
  sim_import_demand(i,r)
;


  sim_price_index(i) ..
    v_sim_W(i) =e= sum(r, v_delta(i,r) ** sigma(i) * price(i,r) ** (1-sigma(i)))  ** (1/(1-sigma(i)));

  sim_import_demand(i,r)  ..
    v_sim_x(i,r) =e= v_utility(i) * v_delta(i,r) ** sigma(i) * (price(i,r)/v_sim_W(i)) ** (-sigma(i))
      + v_mu(i,r);



model CES_sim /sim_price_index, sim_import_demand/;


*  fix calibrated parameters

  v_delta.fx(i,r) = delta(i,r,"modified");
  v_mu.fx(i,r)    = mu(i,r,"modified");


*
*   --- test 1: test the calibrated CES (w. commitment term) in the expected point
*

*  initialize variables
  v_sim_W.l(i)   = 1;
  v_sim_x.l(i,r) = 1;

  v_sim_W.lo(i)  = 0.1;


*  price assumption
  price(i,r)     = p_calib(i,r);

CES_sim.solprint = 1;
solve CES_sim using CNS;
if(CES_sim.numinfes ne 0, abort "problem with the test simulation model");

* -- test if sthe solution matches the expectation
if(sum((i), abs(v_sim_x.l(i,"r") - x_calib(i,"r"))) gt 1.E-2,
  abort "problem with reproducing the assumed calibration point ", v_sim_x.l, x_calib;
 );
*$exit



*
*   --- Note that the import quantity from region "n"
*       will not be the same as x_calib!
*       The reasoning is the following:
*       - fixing delta.fx = 1 for this region implies that
*         the price index will differ from the originally set 1
*       - under the assumption of fix utility this implies
*         a different total expenditure (Y)
*       - although the value share are reproduced, the values themselves not
*       - dividing the import values with the fix price we get
*         a different quantity compared to x_calib
*



*
*   ---  check value shares and total expenditure in the test run
*        The value shares are not fully reproduced, as the
*        calibration is made exact to the x_calib quantities.
*
parameter
           p_check_vs(i,r)  "simulated value shares"
           p_check_Y        "check simulated total expenditure"
;
 p_check_Y       = sum((i,r), v_sim_x.l(i,r) * price(i,r));
 p_check_vs(i,r) = v_sim_x.l(i,r) * price(i,r) / p_check_Y;

*abort "check", p_check_Y, p_check_vs, Y, share, v_sim_x.l, x_calib;




*
*   --- test 2: test the calibrated CES (w. commitment term) in the observed point
*

*  initialize variables
  v_sim_W.l(i)   = 1;
  v_sim_x.l(i,r) = 1;

*  price assumption
  price(i,r)     = p_null(i,r);

CES_sim.solprint = 1;
solve CES_sim using CNS;
if(CES_sim.numinfes ne 0, abort "problem with the test simulation model");

* -- test if sthe solution matches the expectation
if(sum((i), abs(v_sim_x.l(i,"r") - x_null(i,"r"))) gt 1.E-2,
  abort "problem with reproducing the observed calibration point ", v_sim_x.l, x_calib;
 );
*$exit



*
*   ---  Part II.
*
*   ---  Sensitivity analysis with the scenario solver
*        The relative price of imports from region "r"
*        are varied in a range.

*  number of price experiments
$setlocal N 100

set scen "scenarios" /s1*s%N%/;

parameter     p_results(scen, i, r, *, *)   "reporing parameter";


parameter
  ps_price(scen, i, r)    "prices in the experiments"
  ps_x(scen, i, r)        "import demands"
;

*
*   --- Relative price of imports from region "r" are set in [0.05,0.5]
*       The setlocal N defines the number of observations
*
  ps_price(scen, i, "r") = .45 * 1/(%N%-1) * (ord(scen)-1) + .01;
  ps_price(scen, i, "n") = 1;

set scen_dict   "scenario dictionary (for the GUSS solver option)"
/
 scen    .   scenario .     ''
 price   .   param    .     ps_price
 v_sim_x .   level    .     ps_x
/;




*  standard Armington approach

  v_delta.fx(i,r) = delta(i,r,"standard");
  v_mu.fx(i,r)    = mu(i,r,"standard");
  v_sim_W.l(i)    = 1;
  v_sim_x.l(i,r)  = 1;

solve CES_sim using CNS scenario scen_dict;


p_results(scen, i, r, "p", "standard") = ps_price(scen,i,r);
p_results(scen, i, r, "x", "standard") = ps_x(scen,i,r);

  option p_results:4:3:2;
*  abort "check scenario results for the standard case: ", p_results;


*  Modified Armington approach

  v_delta.fx(i,r) = delta(i,r,"modified");
  v_mu.fx(i,r)    = mu(i,r,"modified");
  v_sim_W.l(i)   = 1;
  v_sim_x.l(i,r) = 1;

solve CES_sim using CNS scenario scen_dict;


p_results(scen, i, r, "p", "modified") = ps_price(scen,i,r);
p_results(scen, i, r, "x", "modified") = ps_x(scen,i,r);

execute_unload "witzke_et_al.gdx", p_results;


*
*   --- Prepare plot
*

*
*   --- Put data points in a .dat file
*
file datafile /plot.dat/;
put datafile;
loop(scen,
   put p_results(scen, "i1", "r", "x", "standard"):10:2;
   put ' ',p_results(scen, "i1", "r", "p", "standard"):10:2;
   put ' ',p_results(scen, "i1", "r", "x", "modified"):10:2;
   put ' ',p_results(scen, "i1", "r", "p", "modified"):10:2;
   put /;
);
putclose;


*
*   --- Prepare GNUPLOT script
*
file pltfile /plot.plt/;
put pltfile;
putclose
   'set xlabel "import demand"'/
   'set ylabel "relative price of good from region r"'/
   'set title "Comparison of demand functions (standard vs. commitment versions)"'/
*   'set key off'/
   'set xrange [0:100]'/
   'set term png font arial 13'/
   'set output "plot.png"'/

   'plot "plot.dat" using 1:2 title "standard" with lines, "plot.dat" using 3:4 title "commitment" with lines'
;


* sets and parameters for gnuplotxyz
$setlocal gnuplot_path 'S:\util\gnuplot\bin\'

* Use Gnuplot to generate picture
execute 'call %gnuplot_path%gnuplot plot.plt';

* Use mspaint(Windows) to open image file
execute 'mspaint plot.png';

