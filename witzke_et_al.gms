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
      i1       1         .5
;

table   x_null(i, r)   "demand in the base year"

              n          r
      i1     100         1
;

table    p_calib(i,r)  "hypothetical prices in the emerging trade flow (assumption)"

               n           r
      i1       1          .25
;


sigma(i) = 5;
rho(i)   = 1 / sigma(i) - 1;
Y        = 100;

share(i, r, "data") = p_null(i,r) * x_null(i,r) / sum( (ii,rr), p_null(ii,rr) * x_null(ii,rr));


variables
  v_utility(i)                 "utility"
  v_x(i,r,points)              "import demand"
  v_delta(i,r)                 "CES share parameters (for calibration)"
  v_W(i,points)                "price index"
  v_mu(i,r)                    "commitment parameters"
;


equations
  price_index(i, points)       "price index"
  import_demand(i, r, points)  "optimal import demand, F.O.C. of expenditure minimization under fix utility"
;

  price_index(i, points) $ v_W.range(i,points) ..
    v_W(i,points) =e= sum(r, v_delta(i,r) ** sigma(i) * p(i,r,points) ** (1-sigma(i)))  ** (1/(1-sigma(i)));

  import_demand(i,r,points) $ (v_delta.range(i,r) and v_x.l(i,r,points)) ..
    v_x(i,r,points) =e= v_utility(i) * v_delta(i,r) ** sigma(i) * (p(i,r,points)/v_W(i,points)) ** (-sigma(i))
      + v_mu(i,r);


*
*   --- without loss of generality, one of the delta's is set to unity
*       (changing this would only cause a different scaling (monotone transformation)
*        in the utility but do not change the utility ordering)
*

 delta(i, "n", approaches) = 1;

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

*  demand
 v_x.fx(i,r,"data")          = x_null(i,r);
 v_x.fx(i,r,"assumption")    = 0;

*  prices
 p(i,r,"data")               = p_null(i,r);
 p(i,r,"assumption")         = 0;
 v_W.L(i,"data")             = 1;
 v_W.fx(i,"assumption")      = 1;

*  CES parameters (initialization prevents numerical problems)
 v_delta.fx(i,"n")           = 1;
 v_delta.l (i,"r")           = .5;
*  No commitments in the standard case
 mu(i, r, "standard")        = 0;
 v_mu.fx(i,r)                = 0;

 CES.holdfixed = 1;
 CES.solprint = 1;
 solve CES using CNS;

 if(CES.numinfes ne 0, abort "problem with standard calibration");


* save share parameters in the standard case
 delta(i,r,"standard") =  v_delta.L(i,r);

*
*   ---  Calibrate the CES funciton, modified version
*        equation (5) in the paper
*

* -- the committment on goods from country 'n' remains zero
  v_mu.fx(i, 'n') = 0;


 parameter x_calib(i,r) "demand in the second calibration point (assumption)";
 x_calib(i,r) = share(i,r,"assumption") * Y / p_calib(i,r);


  v_x.fx(i,r,"assumption")          = x_calib(i,r);
  p(i,r,"assumption")               = p_calib(i,r);
  v_W.lo(i,"assumption")            = 0;
  v_W.up(i,"assumption")            = +inf;
  v_W.l(i,"assumption")             = 1;


*  allow for non-zero commitment parameters for country 'r'
  v_mu.l(i, 'r')    = 1;
  v_mu.lo(i, 'r')   = -inf;
  v_mu.up(i, 'r')   = +inf;


 CES.solprint = 1;
 solve CES using CNS;
 if(CES.numinfes ne 0, abort "problem with modified calibration");


 delta(i,r,"modified") =  v_delta.L(i,r);
 mu(i,r,"modified")    =  v_mu.L(i,r);

display "check calibration parameters", delta, mu;


*
*   --- Run test scenario
*



*
*   ---   Simulation model
*

parameters
  price(i,r)  "current price in simulation (exogenous)"
;




*
*   ---   Define a simulation model
*         The reason why we need a new model is that the composite price & quantity
*         variables above were dependent on the calibration points ("points").
*
variables
  v_sim_W(i)       "composite price index in simulations"
  v_sim_x(i,r)     "import demand in simulations"
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

*  initialize variables
  v_sim_W.l(i)   = 1;
  v_sim_x.l(i,r) = 1;

*  price assumption
  price(i,"n")     = 1;
  price(i,"r")     = .25;


solve CES_sim using CNS;

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
  ps_price(scen, i, "r") = .45 * 1/(%N%-1) * (ord(scen)-1) + .05;
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



*  Modified Armington approach

  v_delta.fx(i,r) = delta(i,r,"modified");
  v_mu.fx(i,r)    = mu(i,r,"modified");
  v_sim_W.l(i)   = 1;
  v_sim_x.l(i,r) = 1;

solve CES_sim using CNS scenario scen_dict;


p_results(scen, i, r, "p", "modified") = ps_price(scen,i,r);
p_results(scen, i, r, "x", "modified") = ps_x(scen,i,r);

 option p_results:4:3:2;

display p_results;

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
   'set key off'/
   'set term png font arial 13'/
   'set output "plot.png"'/

   'plot "plot.dat" using 1:2 with lines, "plot.dat" using 3:4 with lines'
;


* sets and parameters for gnuplotxyz
$setlocal gnuplot_path 'S:\util\gnuplot\bin\'

* Use Gnuplot to generate picture
execute 'call %gnuplot_path%gnuplot plot.plt';

* Use mspaint(Windows) to open image file
execute 'mspaint plot.png';

