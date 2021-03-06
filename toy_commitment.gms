$ontext
  test different CES functional forms

  1) NOT including productivity multiplier term (gamma) -- version 1 below
  2) including productivity multiplier term (gamma)     -- version 2 below
  3) using a modified CES function form with commitment-terms

  - by including gamma we can force the adding up condition on the deltas
  - in all cases the benchmark utility is scaled to total demand quantity
  - the price index in these cases is scaled to the average price

  References
  ----------

    (standard):  Rutherford, Thomas F. 1995. "CES Demand Functions: Hints and Formulae."

    (commitment): Witzke, Peter, Marcel Adenäuer, Wolfgang Britz, and Thomas Heckelei. 2005.
                  "Modelling EU Sugar Market Scenarios with a Modified Armington Appoach." In IATRC Symposium. Sevilla, Spain.
$offtext

$offlisting

*1) with our without the gamma term (productivity mulitplier)

* two goods case
 set i  "goods" /1*2/;


* input data
 option seed=12234;

 parameter
          pnull(i)           "prices"
          xnull(i)           "demand"
          pindex_null        "price index"
 ;
 pnull(i)    = uniform(0.5,1);
 xnull(i)    = uniform(0.5,1);

 display "check initial price/quantity framework", pnull, xnull;

* set a small demand share for good2
* only relevant for the commitment term-approach
 xnull("2") = 1.E-2;

*
*   --- Put observed quantities into data files for further plotting 
*
file horizontal_obs /horizontal_obs.dat/;
put horizontal_obs;
           put '0';
           put ' ',xnull("2"):10:2;
           put /;
           put '2';
           put ' ',xnull("2"):10:2;
           put /;
putclose;

file vertical_obs /vertical_obs.dat/;
put vertical_obs;
           put xnull("1"):10:2;
           put ' ','0';
           put /;
           put xnull("1"):10:2;
           put ' ','2';
           put /;
putclose;



 parameter
           rho               "substit. param."
           sigma             "subst. elasticity"
           p(i)              "observed prices"

           calib_results(*,*,*)  "calibration results -- reporting"
 ;

 variable
          gamma              "productivity term"
          delta(i)           "share params."
          mu(i)              "commitment terms"
          x(i)               "demand quantities"
          u                  "utility"
          pindex             "price index"
 ;

 equation
          ces1               "utility aggregator, version 1"
          ces2               "utility aggregator, version 2"
          addup_deltas       "adding up condition for deltas"


* dual approach yielding Hicksian demand curves
          demand1(i)         "Hicksian demand curves, version 1"
          demand2(i)         "Hicksian demand curves, version 2"
          index1             "price index, version 1"
          index2             "price index, version 2"

 ;

 ces1 ..   u =e= sum(i, delta(i) * (x(i) - mu(i))**rho)**(1/rho);

 ces2 ..   u =e= gamma * sum(i, delta(i) * (x(i) - mu(i))**rho)**(1/rho);

 addup_deltas .. 1 =e= sum(i, delta(i));

 demand1(i) $ [x.range(i) or delta.range(i)] .. x(i) - mu(i) =e= delta(i)**sigma * (pindex/p(i))** sigma * u;

 demand2(i) .. x(i) - mu(i) =e= (1/gamma) * (delta(i)*gamma)**sigma * (pindex/p(i))** sigma * u;

 index1 ..     pindex =e= sum(i, delta(i)**sigma * p(i)**(1-sigma))**(1/(1-sigma));

 index2 ..     pindex =e= (1/gamma) * sum(i, delta(i)**sigma * p(i)**(1-sigma))**(1/(1-sigma));


*
* --- version 1a)  without gamma term in the CES function
*


* first scale the utility to total consumption quantity in benchmark
* then the price index will be the average price
 model calib1 /demand1, index1/;

 sigma = 4; 
 rho   = (sigma-1)/sigma;

 mu.fx(i)    = 0;

 p(i)    = pnull(i);
 x.fx(i) = xnull(i);

 pindex_null = sum(i, pnull(i) * xnull(i)) / sum(i, xnull(i));
 pindex.l    = pindex_null;
 pindex.lo   = eps;
 u.fx        = sum(i, xnull(i));

 display pindex.l, u.l;

*  initialize
 delta.l (i)   = 1;
 delta.lo(i)   = eps;


 calib1.solprint  = 1;
 solve calib1 using CNS;
* automatic calibration test
 if(abs(pindex.l - pindex_null) gt 1.E-4, abort "price index not recovered");

  calib_results("gamma", " ", "version1")    = 0;
  calib_results("delta", i  , "version1")    = delta.l(i);
  calib_results("utility"," ","version1")    = u.l;
  calib_results("pindex", " ","version1")    = pindex.l;

* check if the original CES utility aggregator reproduces the same value
 model cesagg1/ces1/;

  u.lo        = 0;
  u.up        = +inf;
  x.fx(i)     = x.l(i);
  delta.fx(i) = delta.l(i);


  cesagg1.solprint   = 1;
  cesagg1.iterlim    = 0;
 solve cesagg1 using CNS;
 if(cesagg1.suminfes gt 1.E-4, abort "problem with recovering utility with the CES aggregator");


*
* --- version 2): calibrate a CES aggregator with the gamma term
*                 and enforce the adding-up condition for deltas
*
 model calib2 /demand2, index2, addup_deltas/;


  u.fx          = u.l;
  delta.lo(i)   = eps;
  delta.up(i)   = +inf;
  gamma.l       = 1;

 calib2.solprint   = 1;
 solve calib2 using CNS;

* automatic calibration test
 if(abs(pindex.l - pindex_null) gt 1.E-4, abort "price index not recovered");

   calib_results("gamma", " ", "version2") = gamma.l;
   calib_results("delta", i  , "version2") = delta.l(i);
   calib_results("utility"," ","version2") = u.l;
   calib_results("pindex", " ","version2") = pindex.l;



* check if the original CES utility aggregator reproduces the same value
 model cesagg2 /ces2/;

  u.lo        = 0;
  u.up        = +inf;
  x.fx(i)     = x.l(i);
  delta.fx(i) = delta.l(i);
  mu.fx(i)    = 0;
  gamma.fx    = gamma.l;

  cesagg2.solprint   = 1;
  cesagg2.iterlim    = 0;
 solve cesagg2 using CNS;
 if(cesagg2.suminfes gt 1.E-4, abort "problem with recovering utility with the CES aggregator");



* another variant when one of the deltas is fixed to unity
* then we can not apply the adding up condition


 model calib3 /demand2, index2/;


  u.fx          = u.l;
  delta.lo(i)   = eps;
  delta.up(i)   = +inf;
  delta.fx("1") = 1;
  gamma.l       = 1;
  gamma.lo      = eps;
  gamma.up      = +inf;

 calib3.solprint   = 1;
 solve calib3 using CNS;

* automatic calibration test
 if(abs(pindex.l - pindex_null) gt 1.E-4, abort "price index not recovered");

   calib_results("gamma", " ", "version2mod") = gamma.l;
   calib_results("delta", i  , "version2mod") = delta.l(i);
   calib_results("utility"," ","version2mod") = u.l;
   calib_results("pindex", " ","version2mod") = pindex.l;


* check if the original CES utility aggregator reproduces the same value
  u.lo        = 0;
  u.up        = +inf;
  x.fx(i)     = x.l(i);
  delta.fx(i) = delta.l(i);
  mu.fx(i)    = 0;
  gamma.fx    = gamma.l;

  cesagg2.solprint   = 1;
  cesagg2.iterlim    = 0;
 solve cesagg2 using CNS;
 if(cesagg2.suminfes gt 1.E-4, abort "problem with recovering utility with the CES aggregator");


  display calib_results;

*
* part c) Using a modified CES function form with commitment-terms
*         The approach allows to deal with the small share problem
*         (and with the problem of zero benchmark demand)
*         Demand reactions of initially small demand shares are therefore magnified
*         compared to the standard CES form approach
*
*    -  Note that the above ces1 and ces2 equations already contain a commitment term "mu"
*    -  Calibrating the commitment term requires an additional degree of freedom. This will
*       be provided by the extra information on expected demand under different relative prices


*  test scenario to get a 'standard CES' reaction
*  (Later we will compare the reaction with the commitment CES to this benchmark)
*
*  - relative price of good decreases significantly
*  - here we also compare the reaction with CES version 1 and version 2

 p("2")       = p("2") * .5;
 x.lo(i)      = eps;
 x.up(i)      = +inf;
 pindex.lo    = eps;
 pindex.up    = +inf;
 u.fx         = u.l;

 parameter results(i,*,*)  "simulated results";

 delta.fx(i)  = calib_results("delta", i,   "version1");
 gamma.fx     = calib_results("gamma", " ", "version1");

 model ces_sim1 "simulation model for the first CES approach (without gamma term)" /demand1, index1/;

 ces_sim1.solprint   = 1;
 solve ces_sim1 using CNS;
 if(ces_sim1.numinfes ne 0, abort "problem with simulation CES version 1");

 results(i, "x", "version1") = x.l(i);


* set back demand levels to a default
 x.l(i)       = 1;

 delta.fx(i)  = calib_results("delta", i,"version2");
 gamma.fx     = calib_results("gamma", " ", "version2");

 model ces_sim2 "simulation model for the second CES approach (with gamma term)" /demand2, index2/;

 ces_sim2.solprint   = 1;
 solve ces_sim2 using CNS;
 if(ces_sim2.numinfes ne 0, abort "problem with simulation CES version 2");

 results(i, "x", "version2") = x.l(i);



* set back demand levels to a default
 x.l(i)       = 1;

 delta.fx(i)  = calib_results("delta", i,   "version2mod");
 gamma.fx     = calib_results("gamma", " ", "version2mod");


 ces_sim2.solprint   = 1;
 solve ces_sim2 using CNS;
 if(ces_sim2.numinfes ne 0, abort "problem with simulation CES version 2");

 results(i, "x", "version2mod") = x.l(i);
 display "check simulated results (should be identical with the versions so far)", results;



*
*   ---   Commitment version, see Witzke et al. (2005)
*

 set points "calibration points"  / observed, expected/;

 positive variables
          pcom(i,points)          "prices -- commitment version"
          xcom(i,points)          "import demand -- commitment version"
          pindexcom(points)       "price indexes -- commitment version"
          ucom(points)            "utility -- commitment versions"
 ;

* we include non-zero commitment term for one commodity only
  mu.fx("1")    = 0;
  mu.lo("2")    = -inf;
  mu.up("2")    = 0;
  mu.l ("2")    = -0.5 * xnull("2");

* assume significantly more demand at the same (lowered) relative price for good 2
 pcom.fx(i, "expected")         = p(i);
 pcom.fx(i, "observed")         = pnull(i);
 xcom.fx(i, "expected")         = x.l(i);
 xcom.fx("2", "expected")       = x.l("2") * 5;
 xcom.fx(i, "observed")         = xnull(i);

* free the quantity variables for good 1 (the one with zero commitment term)
 xcom.lo(i, "expected") $ (not mu.range(i))   = eps;
 xcom.up(i, "expected") $ (not mu.range(i))   = +inf;

 option pcom:3:1:1;
 option xcom:3:1:1;
 display "price/quantity setting", pcom.l, xcom.l;

 display xnull;

 equation
         demand_com(i, points)     "import demand equation"
         index_com(points)         "price index"
         ces_com(points)            "utility aggregator"
 ;

 ces_com(points) $ ucom.range(points) ..

     ucom(points) =e= sum(i, delta(i) * (xcom(i,points) - mu(i))**rho)**(1/rho);

 demand_com(i,points) $ [xcom.range(i,points) or delta.range(i)] ..

     xcom(i,points) - mu(i) =e= delta(i)**sigma * (pindexcom(points)/pcom(i,points))** sigma * ucom(points);

 index_com(points) $ pindexcom.range(points)..

     pindexcom(points) =e= sum(i, delta(i)**sigma * pcom(i,points)**(1-sigma))**(1/(1-sigma));


 model calib_commit "calibration model for the CES commitment version" /index_com, demand_com/;

*
*   --- further initialization
*
 delta.lo(i)  = 0;
 delta.up(i)  = +inf;
 delta.l (i)  = calib_results("delta", i,   "version1");

 pindexcom.lo (points)     = 1.E-3;
 pindexcom.l  (points)     = calib_results("pindex", " ",   "version1");

 ucom.fx (points)          = u.l;

 calib_commit.solprint   = 1;
* calib_commit.iterlim    = 0;
 calib_commit.holdfixed  = 0;
 solve calib_commit using CNS;
 if(calib_commit.numinfes ne 0, abort "problem with the calibration of the commitment version");

*
*   --- Test if we are still at the same isoquant
*
  ucom.lo(points)        = 0;
  ucom.up(points)        = +inf;
  xcom.fx(i,points)      = xcom.l(i,points);
  delta.fx(i)            = delta.l(i);
  mu.fx(i)               = mu.l(i);

 model testu_commit "test utility in the commitment version" /ces_com/;

 solve testu_commit using cns;
 if(testu_commit.numinfes ne 0, abort "problem with utility test in the commitment version");
 if(abs(ucom.l("observed") - ucom.l("expected")) gt 1.E-4, abort "not on the same isoquant");




*
*   ---  Run a test simulation with the calibrated commitment version
*

 ucom.fx(points)         = ucom.l(points);

* -- if both the deltas and the import demands are fixed then the import demand equation is swithced off
*    (see the dollar conditionals in the equation definition)
 xcom.lo(i,"expected")  = eps;
 xcom.up(i,"expected")  = +inf;
* -- similarly, we switch of the 'observed' price index equation by fixing the variable
 pindexcom.fx("observed") = pindexcom.l("observed");


solve calib_commit using CNS;
if(calib_commit.numinfes ne 0, abort "problem with the test simulation with the commitment version");

 results(i, "x", "commitment") = xcom.l(i, "expected");

 display "check simulation results in all versions", results;

*
*   --- Put expected quantities into data files for further plotting 
*
file horizontal /horizontal_commit.dat/;
put horizontal;
           put '0';
           put ' ',xcom.l("2","expected"):10:2;
           put /;
           put '2';
           put ' ',xcom.l("2","expected"):10:2;
           put /;
putclose;

file vertical /vertical_commit.dat/;
put vertical;
           put xcom.l("1","expected"):10:2;
           put ' ','0';
           put /;
           put xcom.l("1","expected"):10:2;
           put ' ','2';
           put /;
putclose;

execute_unload "xcommit.gdx", xcom;

*
*   ---  Sensitivity Analysis with the different relative prices
*        -------------------------------------------------------
*
*        The objective is to draw the indifference curve
*        with the commitment version of the demand system.
*

*        number of price experiments
$setlocal N 100

set scen "scenarios for the SCENARIO solver" /s1*s%N%/;

parameter
  ps_price(scen, i, points)    "prices in the experiments -- assumptions"
  ps_x(scen, i, points)        "import demands            -- model results"
;


Set ma "GUSS Model Attributes" / modelstat, solvestat, objval /;
$eval myNomatchLimit %N%*2

Parameter
          o "additional GUSS solver options"     / NoMatchLimit  %myNomatchLimit% /
          r_s(scen,ma)                           "Solution status report -- generated by GUSS"
;

*
*   --- Relative price of imports from region "2" is varied in a range
*
  display  pcom.l;

  ps_price(scen, i, points)       = pcom.l(i, points);
  ps_price(scen, "2", "expected") = pcom.l("2","expected") * [3 * 1/(%N%-1) * (ord(scen)-1)] + 1.E-1;

  option ps_price:3:2:1;
  display ps_price;

set scen_dict   "scenario dictionary (for the GUSS solver option)"
/
 scen    .   scenario .     ''
 pcom    .   fixed    .     ps_price
 xcom    .   level    .     ps_x
 o       .   opt      .     r_s
/;


 solve calib_commit using CNS scenario scen_dict;

*
*   ---  Some results and model statistics
*
 display "Model statistics: ", r_s;


 option ps_x:3:2:1;
* display "solutions for import demand: ", ps_x;

*
*   ---  Prepare one big reporting parameter
*
parameter     p_results(scen, i, points, *)   "reporing parameter";

 p_results(scen, i, points, "price")     = ps_price(scen, i, points);
 p_results(scen, i, points, "demand")    = ps_x(scen, i, points);


 display p_results;


*
*   --- Prepare plot
*

*
*   --- Put data points in a .dat file
*       Note that the 'expected' dimension contains the simulated results
file datafile /plot_commitment.dat/;
put datafile;
loop(scen,
           put p_results(scen, "1", "expected", "demand"):10:2;
           put ' ',p_results(scen, "2", "expected", "demand"):10:2;
           put /;
);
putclose;


*
*   --- Prepare GNUPLOT script
*
file pltfile /plot.plt/;
put pltfile;
putclose
   'set xlabel "import demand from country 1"'/
   'set ylabel "import demand from country 2"'/
   'set title  "Price SA with the commitment version"'/
*   'set key off'/
   'set xrange [0:2]'/
   'set yrange [0:2]'/
   'set term png font arial 13'/
   'set output "plot.png"'/

   'plot "plot_commitment.dat" using 1:2 title "indifference curve" with lines, \' /
   '"horizontal_commit.dat" using 1:2 title "x2 exp." with lines, "vertical_commit.dat" using 1:2 title "x1 exp." with lines, \' /
   '"horizontal_obs.dat" using 1:2 title "x2 obs." with lines, "vertical_obs.dat" using 1:2 title "x1 obs." with lines' 
;


* sets and parameters for gnuplotxyz
$setlocal gnuplot_path 'S:\util\gnuplot\bin\'

* Use Gnuplot to generate picture
execute 'call %gnuplot_path%gnuplot plot.plt';

* Use mspaint(Windows) to open image file
execute 'mspaint plot.png';




