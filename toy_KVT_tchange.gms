$ontext
  modified Kuiper-Tongeren approach
  Basic idea: Hicks non-neutral parameters are a function of relative price changes

  1) NOT including productivity multiplier term (gamma) -- version 1 below
  2) including productivity multiplier term (gamma)     -- version 2 below
  3) using a modified CES function form with a multiplicative term

  - by including gamma we can force the adding up condition on the deltas
  - in all cases the benchmark utility is scaled to total demand quantity
  - the price index in these cases is scaled to the average price

  References
  ----------

    (standard):  Rutherford, Thomas F. 1995. CES Demand Functions: Hints and Formulae.

    (multiplicative calib. term): Kuiper, Marijke, and Frank van Tongeren. 2006.
                                 "Using Gravity to Move Armington -  an Empirical Approach to the Small Initial Trade Share Problem in General Equilibrium Models."
                                 In Presented at the 9th Annual Conference on Global Economic Analysis, Addis Ababa, Ethiopia.

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


* set a small demand share for good2
* only relevant for the commitment term-approach
 xnull("2") = 0.01;

 display "check initial price/quantity framework", pnull, xnull;

 parameter
           rho               "substit. param."
           sigma             "subst. elasticity"
           p(i)              "observed prices"

*  modified KVT approach
          a(i)               "parameter of the linear function defining the change in delta"
          b(i)               "parameter of the linear function defining the change in delta"


           calib_results(*,*,*,*)  "calibration results -- reporting"
 ;

 variable
          gamma              "productivity term"
          delta(i)           "share params."
          x(i)               "demand quantities"
          u                  "utility"
          pindex             "price index"

*  calibration parameters related to the factor-augmented CES form
          Hicks              "Hicks neutral technical change parameters"
          theta(i)           "factor augmenting tech. change (equivalent term with respect to utility function?)"


 ;

 equation
          ces1               "utility aggregator, simple version without gamma term"
          ces2               "utility aggregator, simple version with    gamma term"
          ces_aug            "CES aggregator with factor-augmenting techn. change"
          addup_deltas       "adding up condition for deltas"


* dual approach yielding Hicksian demand curves
          demand1(i)         "Hicksian demand curves, version 1"
          demand2(i)         "Hicksian demand curves, version 2"
          demand_aug(i)      "Hicksian demand, factor-aug. version"
          index1             "price index, version 1"
          index2             "price index, version 2"
          index_aug          "price index, factor-aug. version"

          Hicks_neutral      "Hicks neutral technological change"

*  modified KVT
          transition(i)       "transition function driving the preference-change parameters"

 ;

 ces1 ..      u =e= sum(i, delta(i) * x(i)**rho)**(1/rho);

 ces2 ..      u =e= gamma * sum(i, delta(i) * x(i)**rho)**(1/rho);

 ces_aug ..   u =e= (gamma * Hicks) * sum(i, delta(i) * ((x(i)*theta(i))**rho))**(1/rho);

 addup_deltas .. 1 =e= sum(i, delta(i));

 demand1(i) $ [x.range(i) or delta.range(i)] .. x(i) =e= delta(i)**sigma * (pindex/p(i))** sigma * u;



*
*   ---     Note that as we included the gamma term in the price index
*           we have to correct the demand equation with gamma**(sigma-1)
*
 demand2(i) ..         x(i) =e= (gamma**(sigma-1)) * (delta(i)**sigma) * (pindex/p(i))** sigma * u;

 demand_aug(i) ..      x(i) =e= (gamma**(sigma-1)) * (delta(i)**sigma) * (theta(i)**(rho*sigma)) * (pindex/p(i))**sigma * u;

 index1 ..           pindex =e= sum(i, delta(i)**sigma * p(i)**(1-sigma))**(1/(1-sigma));

 index2 ..           pindex =e= (1/gamma) * sum(i, delta(i)**sigma * p(i)**(1-sigma))**(1/(1-sigma));

 index_aug ..        pindex =e= (1/gamma) * sum(i, (delta(i)**sigma) * (theta(i)**(sigma-1)) * (p(i)**(1-sigma)))**(1/(1-sigma));

 Hicks_neutral ..     gamma =e= Hicks * gamma.l;

 transition(i) $ theta.range(i)  ..  theta(i) =e= b(i) * [(p(i) / pindex) / ( pnull(i)/pindex_null)]
                                               + a(i) * [(p(i) / pindex) / ( pnull(i)/pindex_null) - 1] ;

*
* --- version 1a)  without gamma term in the CES function
*


* first scale the utility to total consumption quantity in benchmark
* then the price index will be the average price
 model calib1 /demand1, index1/;

 sigma = 4;
 rho   = (sigma-1)/sigma;

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
 if(abs(pindex.l - pindex_null) gt 1.E-4, abort "price index not recovered", pindex_null);

  calib_results("gamma", " ", "version1","observed")    = 0;
  calib_results("delta", i  , "version1","observed")    = delta.l(i);
  calib_results("utility"," ","version1","observed")    = u.l;
  calib_results("pindex", " ","version1","observed")    = pindex.l;

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
 if(abs(pindex.l - pindex_null) gt 1.E-4, abort "price index not recovered", pindex_null);

   calib_results("gamma", " ", "version2","observed") = gamma.l;
   calib_results("delta", i  , "version2","observed") = delta.l(i);
   calib_results("utility"," ","version2","observed") = u.l;
   calib_results("pindex", " ","version2","observed") = pindex.l;



* check if the original CES utility aggregator reproduces the same value
 model cesagg2 /ces2/;

  u.lo        = 0;
  u.up        = +inf;
  x.fx(i)     = x.l(i);
  delta.fx(i) = delta.l(i);
  gamma.fx    = gamma.l;

  cesagg2.solprint   = 1;
  cesagg2.iterlim    = 0;
 solve cesagg2 using CNS;
 if(cesagg2.suminfes gt 1.E-4, abort "problem with recovering utility with the CES aggregator");


*
* --- version 2 mod): another variant when one of the deltas is fixed to unity
*                     then we can not apply the adding up condition
*
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
 if(calib3.numinfes ne 0, abort "problem with calibrating version 2mod");

* automatic calibration test
 if(abs(pindex.l - pindex_null) gt 1.E-4, abort "price index not recovered");

   calib_results("gamma", " ", "version2mod","observed") = gamma.l;
   calib_results("delta", i  , "version2mod","observed") = delta.l(i);
   calib_results("utility"," ","version2mod","observed") = u.l;
   calib_results("pindex", " ","version2mod","observed") = pindex.l;

*
*   --- test if the original CES utility aggregator reproduces the same value
*
  u.lo        = 0;
  u.up        = +inf;
  x.fx(i)     = x.l(i);
  delta.fx(i) = delta.l(i);
  gamma.fx    = gamma.l;

  cesagg2.solprint   = 1;
  cesagg2.iterlim    = 0;
 solve cesagg2 using CNS;
 if(cesagg2.suminfes gt 1.E-4, abort "problem with recovering utility with the CES aggregator");


*
*   --- Calibration of version 3)
*                  CES form with Hicks-neutral and augmented technological change
*                  Here we set the technological change parameters to unity
*                  i.e. no technological change is assumed so far
*
 model ces_tchange / index_aug, demand_aug/;

  Hicks.fx      =  1;
  theta.fx(i)   =  1;
  u.fx          =  u.l;

  delta.l(i)    = calib_results("delta", i,   "version2mod","observed");
  delta.lo(i)   = eps;
  delta.up(i)   = +inf;
  delta.fx("1") = 1;

  gamma.l       = calib_results("gamma", " ", "version2mod","observed");
  gamma.lo      =  0;
  gamma.up      =  +inf;

 ces_tchange.solprint   = 1;
 solve ces_tchange using CNS;
 if(ces_tchange.numinfes ne 0, abort "problem with simulation CES with technological change");

   calib_results("gamma", " ", "tchange","observed") = gamma.l;
   calib_results("delta", i  , "tchange","observed") = delta.l(i);
   calib_results("utility"," ","tchange","observed") = u.l;
   calib_results("pindex", " ","tchange","observed") = pindex.l;
   calib_results("theta", i  ,"tchange","observed")  = theta.l(i);
   calib_results("Hicks", " ", "tchange","observed") = Hicks.l;

*
*   --- test if the original CES utility aggregator reproduces the same value
*
 model test_cesaug /ces_aug/;

  u.lo        = 0;
  u.up        = +inf;
  x.fx(i)     = x.l(i);
  delta.fx(i) = delta.l(i);
  gamma.fx    = gamma.l;

  test_cesaug.solprint   = 1;
  test_cesaug.iterlim    = 0;
 solve test_cesaug using CNS;
 if(test_cesaug.suminfes gt 1.E-4, abort "problem with recovering utility with the CES aggregator");


 option calib_results:3:3:1;
 display calib_results;

*
*  --- test scenario to get 'standard CES' reactions
*      (versions 1 - 3)
*      Later we will compare the reaction with the commitment CES to this benchmark
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



*
*   ---   Test simulation version 1
*
 delta.fx(i)  = calib_results("delta", i,   "version1","observed");
 gamma.fx     = calib_results("gamma", " ", "version1","observed");

 model ces_sim1 "simulation model for the first CES approach (without gamma term)" /demand1, index1/;

 ces_sim1.solprint   = 1;
 solve ces_sim1 using CNS;
 if(ces_sim1.numinfes ne 0, abort "problem with simulation CES version 1");

 results(i, "x", "version1") = x.l(i);


*
*   ---   Test simulation version 2
*
* set back demand levels to a default
 x.l(i)       = 1;

 delta.fx(i)  = calib_results("delta", i,"version2","observed");
 gamma.fx     = calib_results("gamma", " ", "version2","observed");

 model ces_sim2 "simulation model for the second CES approach (with gamma term)" /demand2, index2/;

 ces_sim2.solprint   = 1;
 solve ces_sim2 using CNS;
 if(ces_sim2.numinfes ne 0, abort "problem with simulation CES version 2");

 results(i, "x", "version2") = x.l(i);


*
*   ---   Test simulation version 2 mod
*
* set back demand levels to a default
 x.l(i)       = 1;

 delta.fx(i)  = calib_results("delta", i,   "version2mod","observed");
 gamma.fx     = calib_results("gamma", " ", "version2mod","observed");


 ces_sim2.solprint   = 1;
 solve ces_sim2 using CNS;
 if(ces_sim2.numinfes ne 0, abort "problem with simulation CES version 2");

 results(i, "x", "version2mod") = x.l(i);


*
*   ---   Test simulation version 3 (techn. change)
*
 model test_tchange / index_aug, demand_aug, Hicks_neutral/;
* set back demand levels to a default
  x.l(i)       = 1;


 Hicks.fx     = 1;
 theta.fx(i)  = 1;
 delta.fx(i)  = calib_results("delta", i,   "tchange","observed");

 gamma.lo     = 0;
 gamma.l      = calib_results("gamma", " ", "tchange","observed");
 gamma.up     = gamma.l * 10;


 test_tchange.solprint   = 1;
 solve test_tchange using CNS;
 if(test_tchange.numinfes ne 0, abort "problem with simulation CES version 2");

 results(i, "x", "tchange") = x.l(i);

 display "check simulated results (should be identical with the versions so far)", results;

*
*   --- Put test simulation result (showing the small share problem)
*       into data files for further plotting
*
file point_small /point_small.dat/;
put point_small;
           put x.l("1"):10:5;
           put ' ',x.l("2"):10:5;
           put /;
putclose;


*$exit
*
*   ---   Kuiper - van Tongeren (2006)
*
*
* part c) Using a modified CES function form with multiplicative calibration terms
*         The approach allows to deal with the small share problem
*         Demand reactions of initially small demand shares are therefore magnified
*         compared to the standard CES form approach
*
*    -  The multiplicative term is a factor-augmented technological change term (theta)
*    -  The CES specification also contains the Hicks-neutral techn. change term, but it is kept at constant 1
*
*    (1) The functional form is first calibrated to the observed point as in version 3 above, with theta.fx = 1
*    (2) Then we fix the deltas and calibrate with the theta's to the expected point


*
*  --- step 2 from above: Calibrate to expected point
*


* Expected point
* ---------------
*
* Assume significantly more demand at the same (lowered) relative price for good 2
* [That's a textbook case for the small share problem, i.e. we would expect
* more reaction than made possible by the standard CES functional form]
*
*
*  --- A difficulty with the Kuiper-Tongeren approach arises:
*      With how much would the demand for imports from country 1 decrease?
*      It is normally defined by gravity estimates. Here we need to provide
*      an estimate that would keep utility at a more-or-less constant level.
*
*      One solution is to use the import demand from country 1 calculated by the commitment version.
*
 parameter   expected(i,*);

 x.fx(i)                   = x.l(i);
 x.fx("2")                 = x.l("2") * 5;


$ontext
*
*  --- use the same expected point as achieved by
*      the commitment version
*
variable xcom;
$gdxin xcommit
$load xcom
$gdxin
 x.fx("1")  =  xcom.l("1","expected");
$offtext

* -- alternatively, simply use a guesstimate for the decrease
 x.fx("1")  =  x.l("1") / 2;

*
*   --- Price development is the same as in the test scenarios above
*       (demand for good 2 is increased at the same relative prices)
*

 display "price/quantity framework in the expected calibration point", p,x.l;

 expected(i,"price")    =  p(i);
 expected(i,"demand")   =  x.l(i);


*
*   ---  Calibrate the factor-augmenting CES specification to the expected point
*        (CES version tchnage, i.e. with Hicks non-neutral theta params)
*        This calculates (1) expected price index and (2) thetas that we need to achive with the transition equation
*
*   ---  Note that the calibration is not independent from the first calibration to the observed point:
*        -  we keep utility at the prevoius level (stay on the same isoquant)

 pindex.l    = sum(i, p(i) * x.l(i)) / sum(i, x.l(i));
 pindex.lo   = eps;
 pindex.up   = 3.0 * pindex.l;

 u.fx          = calib_results("utility"," ","tchange","observed");
 gamma.fx      = calib_results("gamma"," ","tchange","observed");
 delta.fx(i)   = calib_results("delta",i,"tchange","observed");

  Hicks.fx      =  1;

*  initialize
  theta.l (i)   =  1;
  theta.lo(i)   = eps;
  theta.up(i)   = 10.0 * theta.l(i);

 ces_tchange.solprint  = 1;
 solve ces_tchange using CNS;
* automatic calibration test
 if(ces_tchange.suminfes gt 1.E-4, abort "augmented CES could not be calibrated to expected point", calib_results);

   calib_results("gamma", " ", "tchange","expected") = gamma.l;
   calib_results("delta", i  , "tchange","expected") = delta.l(i);
   calib_results("utility"," ","tchange","expected") = u.l;
   calib_results("pindex", " ","tchange","expected") = pindex.l;
   calib_results("theta", i  ,"tchange","expected")  = theta.l(i);
   calib_results("price", i  ,"tchange","expected")  = p(i);
   calib_results("demand", i  ,"tchange","expected") = x.l(i);
   calib_results("Hicks", " ", "tchange","expected") = Hicks.l;

*
*   --- Calibration of the linear price difference function (modified KVT approach)
*

  b(i) = 1;

  a(i) =  [calib_results("theta", i  , "tchange","expected")
           - ((p(i)/pindex.l) / (pnull(i)/pindex_null))]
           /
            [(p(i)/pindex.l) / (pnull(i)/pindex_null) - 1];

 display "Parameters of the transition function: ", a;

*
*
*   --- check if we are at the same isoquant
*
 model isoquant_modKVT / ces_aug, transition/;
  u.lo        = 0;
  u.up        = +inf;
  x.fx(i)     = x.l(i);
  pindex.fx   = pindex.l;
  theta.fx(i) = theta.l(i);

  isoquant_modKVT.solprint   = 1;
  isoquant_modKVT.iterlim    = 0;
 solve isoquant_modKVT using CNS;
 if(isoquant_modKVT.suminfes gt 1.E-4, abort "problem with recovering utility with the CES aggregator");


*
*   ---  Test simulation with the modified KVT approach
*        a) does it reproduce the observed point?

model mod_KVT /index_aug, demand_aug, transition/;


*
*   --- theta is now defined by the transition function
*
 theta.l (i)   = calib_results("theta", i,   "tchange","observed");
 theta.lo(i)   = eps;
 theta.up(i)   = 10 * theta.l(i);

 gamma.fx      = calib_results("gamma", " ", "tchange","observed");

 x.l (i)       = xnull(i);
 x.lo(i)       = eps;
 x.up(i)       = 10 * x.l(i);

 u.fx          = u.l;

 p(i)          = pnull(i);
 pindex.l      = pindex_null;
 pindex.lo     = eps;
 pindex.up     = 10 * pindex.l;

 mod_KVT.solprint   = 1;
 mod_KVT.iterlim    = 0;
 solve mod_KVT using CNS;
 if(mod_KVT.numinfes ne 0, abort "problem with replicating benchmark with the MODIFIED KVT");

*
*   ---  b) does it reproduce the expected point?
*
 theta.l (i)   = calib_results("theta", i,   "tchange","expected");
 theta.lo(i)   = eps;
 theta.up(i)   = 10 * theta.l(i);

 gamma.fx      = calib_results("gamma", " ", "tchange","expected");

 x.l (i)       = calib_results("demand", i,   "tchange","expected");
 x.lo(i)       = eps;
 x.up(i)       = 10 * x.l(i);

 u.fx          = u.l;

 p(i)          = calib_results("price", i,   "tchange","expected");

 pindex.l      = calib_results("pindex"," ", "tchange","expected");
 pindex.lo     = eps;
 pindex.up     = 10 * pindex.l;

 mod_KVT.solprint   = 1;
 mod_KVT.iterlim    = 0;
 solve mod_KVT using CNS;
 if(ces_sim1.numinfes ne 0, abort "problem with replicating benchmark with the MODIFIED KVT");


*
*   --- Put expected quantities into data files for further plotting
*

file point_exp /point_exp.dat/;
put point_exp;
           put x.l("1"):10:5;
           put ' ',x.l("2"):10:5;
           put /;
putclose;

file point_obs /point_obs.dat/;
put point_obs;
           put xnull("1"):10:5;
           put ' ',xnull("2"):10:5;
           put /;
putclose;

*
*   ---  Inspect simulation behavior:
*        Sensitivity Analysis regarding different relative prices
*        -------------------------------------------------------
*
*        The objective is to draw the indifference curve
*        with the modified Kuiper-van Tongeren approach
*


*    fix functional parameters
   gamma.fx        =  gamma.l;
   delta.fx(i)     =  delta.l(i);

*        number of price experiments
$setlocal N 100

set scen "scenarios for the SCENARIO solver" /s1*s%N%/;

parameter
  ps_price(scen, i)    "prices in the experiments -- assumptions"
  ps_x(scen, i)        "import demands            -- model results"
;


Set ma "GUSS Model Attributes" / modelstat, solvestat, objval /;
$eval myNomatchLimit %N%*2

Parameter
          o "additional GUSS solver options"     / NoMatchLimit  %myNomatchLimit% /
          r_s(scen,ma)                           "Solution status report -- generated by GUSS"
;


  ps_price(scen, i)        = p(i);
  ps_price(scen, "2")      = p("2") * [1.5 * 1/(%N%-1) * (ord(scen)-1)] + 2.5 * 1.E-1;

  option ps_price:3:1:1;
  display ps_price;

set scen_dict   "scenario dictionary (for the GUSS solver option)"
/
 scen    .   scenario .     ''
 p       .   param    .     ps_price
 x       .   level    .     ps_x
 o       .   opt      .     r_s
/;

 mod_KVT.iterlim = 1000;

* solve mod_KVT using CNS;
* display calib_results;

*
*   ---  Start the SA with the GUSS solver
*
 solve mod_KVT using CNS scenario scen_dict;

*
*   ---  Some results and model statistics
*
 display "Model statistics: ", r_s;


 option ps_x:3:1:1;
 display "solutions for import demand: ", ps_x;



*
*   ---  Delete results without a modelstat 16 or 15
*
  ps_x(scen, i)     $ (r_s(scen,"modelstat") lt 15) = 0;
  ps_price(scen, i) $ (r_s(scen,"modelstat") lt 15) = 0;

*
*   ---  Prepare one big reporting parameter
*
parameter     p_results(scen, i, *)   "reporing parameter";

 p_results(scen, i,  "price")     = ps_price(scen, i);
 p_results(scen, i,  "demand")    = ps_x(scen, i);


 display p_results;

*
*   --- Prepare plot
*

*
*   --- Put data points in a .dat file
*       Note that the 'expected' dimension contains the simulated results
file datafile /plot_modKVT.dat/;
put datafile;
loop(scen $ p_results(scen,"1","demand"),
           put p_results(scen, "1", "demand"):10:2;
           put ' ',p_results(scen, "2", "demand"):10:2;
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
   'set title  "modified KVT approach"'/
*   'set key off'/
   'set xrange [0:2]'/
   'set yrange [0:2]'/
   'set term png font arial 13'/
   'set output "plot.png"'/
   'set style line 1 lc rgb "black" pt 5'/
   'set style line 2 lc rgb "black" pt 7'/
   'set style line 3 lc rgb "black" pt 9'/
   'set style line 10 linetype 1 lc rgb "black" lw 2'/

   'plot "plot_modKVT.dat" using 1:2 title "indifference curve" with lines ls 10, \' /
   '"point_obs.dat" using 1:2 title "observed point" w p ls 1, \' /
   '"point_exp.dat" using 1:2 title "expected point" w p ls 2, \' /
   '"point_small.dat" using 1:2 title "small-share point" w p ls 3'

;


* sets and parameters for gnuplotxyz
$setlocal gnuplot_path 'S:\util\gnuplot\bin\'

* Use Gnuplot to generate picture
execute 'call %gnuplot_path%gnuplot plot.plt';

* Use mspaint(Windows) to open image file
execute 'mspaint plot.png';


