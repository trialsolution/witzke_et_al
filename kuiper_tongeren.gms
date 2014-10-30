********************************************************************************
$ontext

   GAMS file : KUIPER_TONGEREN.GMS

   @purpose  :
               Objective
               ---------

               Implement a different modified Armington approach using
               the numerical example in Witzke et al. conference paper
               "Modelling EU sugar market scenarios
               with a modified Armington approach".

               This different approach is the one of Kuiper-Tongeren (2006)
               dealing with the small share problem (originally in GTAP)

               One of the issues is to convert their linearized formulas
               back to the level form.


               References
               ----------

               Kuiper, Marijke, and Frank van Tongeren. 2006.
               “Using Gravity to Move Armington -  an Empirical Approach to the Small Initial Trade Share Problem in General Equilibrium Models.”
               In Presented at the 9th Annual Conference on Global Economic Analysis, Addis Ababa, Ethiopia.

               Witzke, Peter, Marcel Adenäuer, Wolfgang Britz, and Thomas Heckelei. 2005.
               “Modelling EU Sugar Market Scenarios with a Modified Armington Appoach.” In IATRC Symposium. Sevilla, Spain.


   @author   :
   @date     : 29.10.14
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
  delta(i, r, approaches)       "CES share parameter (storage)"
  theta(i, r, approaches)       "factor augmenting t. change parameter (storage)"
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
      i1     100         2
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
  v_H(i)                "Hicks neutral technical change parameters"
  v_theta(i,r)                 "factor augmenting tech. change (equivalent term with respect to utility function?)"
  v_bigtheta1(i,r)
  v_bigtheta2(i,r)
;


equations
  price_index(i, points)       "price index"
  import_demand(i, r, points)  "optimal import demand, F.O.C. of expenditure minimization under fix utility"
  big_theta1(i,r)
  big_theta2(i,r)
;

  price_index(i, points) $ v_W.range(i,points) ..
    v_W(i,points) * v_H(i) =e= sum(r, v_bigtheta1(i,r)
                                          * p(i,r,points)**(1-sigma(i)))  ** (1/(1-sigma(i)));

  import_demand(i,r,points) $ (v_delta.range(i,r) and v_x.l(i,r,points)) ..
    v_x(i,r,points) =e= v_H(i)**(sigma(i)-1)
                      * v_utility(i) * v_bigtheta2(i,r)
                      * (p(i,r,points)/v_W(i,points)) ** (-sigma(i));

  big_theta1(i,r) $ (v_theta.range(i,r) or v_delta.range(i,r)) ..
    v_bigtheta1(i,r) =e= (v_delta(i,r)**sigma(i)) * (v_theta(i,r)**(sigma(i) - rho(i)));

  big_theta2(i,r) $ (v_theta.range(i,r) or v_delta.range(i,r)) ..
    v_bigtheta2(i,r) =e= (v_delta(i,r)**sigma(i)) * v_theta(i,r)**(rho(i) * sigma(i));


 model CES /import_demand, price_index, big_theta1, big_theta2 /;



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




*  demand
 v_x.fx(i,r,"data")          = x_null(i,r);
 v_x.fx(i,r,"assumption")    = 0;

*  prices
 p(i,r,"data")               = p_null(i,r);
 p(i,r,"assumption")         = 0;
 v_W.L(i,"data")             = 1;

*  by fixing the price index we switch off the price index equation
*  for the 'assumed' calibration point (not included in the standard case)
 v_W.fx(i,"assumption")      = 1;

*  CES share parameters are non-negative
 v_delta.lo(i,r)             = eps;

*  CES parameters (initialization prevents numerical problems)
 v_delta.fx(i,"n")           = 1;
 v_delta.l (i,"r")           = .5;


*  Hicks neutral technical change parameter fixed to 1 in the standard case
 v_H.fx(i)            = 1;
*  factor augmenting tech. change fixed to 1 in the standard case
 v_theta.fx(i,r)             = 1;

 v_bigtheta1.fx(i,r) = v_delta.l(i,r)**sigma(i) * v_theta.l(i,r)**(sigma(i) - rho(i));
 v_bigtheta2.fx(i,r) = (v_delta.l(i,r) * v_theta.l(i,r)**rho(i)) ** sigma(i);

* free the bigtheta variables for region "r"
 v_bigtheta1.lo(i,"r")   =  0;
 v_bigtheta1.up(i,"r")   =  +inf;
 v_bigtheta2.lo(i,"r")   =  0;
 v_bigtheta2.up(i,"r")   =  +inf;


 CES.holdfixed = 1;
 CES.solprint = 1;
 solve CES using CNS;

 if(CES.numinfes ne 0, abort "problem with standard calibration");
*$exit

* save share parameters in the standard case
 delta(i,r,"standard") =  v_delta.L(i,r);

*
*   ---  Calibration to expected (non-zero) trade with region "r"
*


 parameter x_calib(i,r) "demand in the second calibration point (assumption)";
 x_calib(i,r) = share(i,r,"assumption") * Y / p_calib(i,r);

*  assumed price/quantity point
  v_x.fx(i,r,"assumption")          = x_calib(i,r);
  p(i,r,"assumption")               = p_calib(i,r);

*  free the price index for the 'assumed' calibration point
  v_W.lo(i,"assumption")            = 0;
  v_W.up(i,"assumption")            = +inf;
  v_W.l(i,"assumption")             = 1;


*  the non-neutral techn. change parameter gives the missing degree of freedom
*  for the calibration
  v_theta.lo(i,"r")           = 0.01;
  v_theta.up(i,"r")           = +inf;

 v_delta.l (i,"r")           = .5;

 CES.holdfixed = 0;
* CES.iterlim   = 0;
 CES.solprint  = 1;
 solve CES using CNS;
 if(CES.numinfes ne 0, abort "problem with calibration model (commitment version)");

*$exit

 delta(i,r,"modified") =  v_delta.l(i,r);
 theta(i,r,"modified") =  v_theta.l(i,r);
*display "check calibration parameters", delta, mu, rho, sigma;

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
variables
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
    v_sim_W(i) * v_H(i) =e= sum(r, v_bigtheta1(i,r) * price(i,r) ** (1-sigma(i)))  ** (1/(1-sigma(i)));

  sim_import_demand(i,r)  ..
    v_sim_x(i,r) =e= v_H(i) ** (sigma(i)-1)
                     * v_utility(i) * v_bigtheta2(i,r) * (price(i,r)/v_sim_W(i)) ** (-sigma(i));



model CES_sim /sim_price_index, sim_import_demand, big_theta1, big_theta2/;


*  fix calibrated parameters

  v_delta.fx(i,r)      = delta(i,r,"modified");
  v_theta.fx(i,r)      = v_theta.l(i,r);
  v_bigtheta1.fx(i,r)  = v_bigtheta1.l(i,r);
  v_bigtheta2.fx(i,r)  = v_bigtheta2.l(i,r);



*  initialize variables
  v_sim_W.l(i)   = 1;
  v_sim_x.l(i,r) = 1;

*  price assumption
  price(i,"n")     = 1;
  price(i,"r")     = .25;

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

*abort "check", p_check_Y, p_check_vs, Y, share;




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
  v_theta.fx(i,r) = 1;
  v_sim_W.l(i)    = 1;
  v_sim_x.l(i,r)  = 1;

solve CES_sim using CNS scenario scen_dict;


p_results(scen, i, r, "p", "standard") = ps_price(scen,i,r);
p_results(scen, i, r, "x", "standard") = ps_x(scen,i,r);



*  Modified Armington approach

  v_delta.fx(i,r) = delta(i,r,"modified");
  v_theta.fx(i,r) = theta(i,r,"modified");
  v_sim_W.l(i)   = 1;
  v_sim_x.l(i,r) = 1;

solve CES_sim using CNS scenario scen_dict;


p_results(scen, i, r, "p", "modified") = ps_price(scen,i,r);
p_results(scen, i, r, "x", "modified") = ps_x(scen,i,r);

 option p_results:4:3:2;

display p_results;

execute_unload "kuiper_tongeren.gdx", p_results;




*
*   --- Prepare plot
*

*
*   --- Put data points in a .dat file
*
file datafile /plot_kuiper.dat/;
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
*   --- Prepare GNUPLOT_KUIPER script
*
file pltfile /plot_kuiper.plt/;
put pltfile;
putclose
   'set xlabel "import demand"'/
   'set ylabel "relative price of good from region r"'/
   'set title "Comparison of demand functions (standard vs. commitment versions)"'/
*   'set key off'/
   'set xrange [0:100]'/
   'set term png font arial 13'/
   'set output "plot_kuiper.png"'/

   'plot "plot_kuiper.dat" using 1:2 title "standard" with lines, "plot_kuiper.dat" using 3:4 title "commitment" with lines'
;


* sets and parameters for gnuplot_kuiperxyz
$setlocal gnuplot_path 'D:\util\gnuplot\bin\'

* Use Gnuplot_kuiper to generate picture
execute 'call %gnuplot_path%gnuplot plot_kuiper.plt';

* Use mspaint(Windows) to open image file
execute 'mspaint plot_kuiper.png';

