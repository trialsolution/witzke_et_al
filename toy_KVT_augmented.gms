$ontext
  modified Kuiper-Tongeren approach
  CES version with augmented technological change parameters
  [to be included in an extended form in CAPRI]
  Basic idea: Hicks non-neutral parameters are a function of relative price changes

  3) using a modified CES function form with a multiplicative term

  - by including gamma we can force the adding up condition on the deltas
  - in all cases the benchmark utility is scaled to total demand quantity
  - the price index in these cases is scaled to the average price
  - the original KVT approach is extended with a transition function => that enables benchmark replication

  References
  ----------

    (standard):  Rutherford, Thomas F. 1995. CES Demand Functions: Hints and Formulae.

    (multiplicative calib. term): Kuiper, Marijke, and Frank van Tongeren. 2006.
                                 "Using Gravity to Move Armington -  an Empirical Approach to the Small Initial Trade Share Problem in General Equilibrium Models."
                                 In Presented at the 9th Annual Conference on Global Economic Analysis, Addis Ababa, Ethiopia.

$offtext



*
*
*   --- setting for gnuplot
*
$setlocal gnuplot_path 'C:\Users\himicmi\Downloads\gnuplot\bin\'


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
          ces_aug            "CES aggregator with factor-augmenting techn. change"
          addup_deltas       "adding up condition for deltas"


* dual approach yielding Hicksian demand curves
          demand_aug(i)      "Hicksian demand, factor-aug. version"
          index_aug          "price index, factor-aug. version"

          Hicks_neutral      "Hicks neutral technological change"

*  modified KVT
          transition(i)       "transition function driving the preference-change parameters"

 ;

 ces_aug ..   u =e= (gamma * Hicks) * sum(i, delta(i) * ((x(i)*theta(i))**rho))**(1/rho);

 addup_deltas .. 1 =e= sum(i, delta(i));

*
*   ---     Note that as we included the gamma term in the price index
*           we have to correct the demand equation with gamma**(sigma-1)
*
 demand_aug(i) ..      x(i) =e= (gamma**(sigma-1)) * (delta(i)**sigma) * (theta(i)**(rho*sigma)) * (pindex/p(i))**sigma * u;

 index_aug ..        pindex =e= (1/gamma) * sum(i, (delta(i)**sigma) * (theta(i)**(sigma-1)) * (p(i)**(1-sigma)))**(1/(1-sigma));

 Hicks_neutral ..     gamma =e= Hicks * gamma.l;

 transition(i) $ theta.range(i)  ..  theta(i) =e= b(i) * [(p(i) / pindex) / ( pnull(i)/pindex_null)]
                                               + a(i) * [(p(i) / pindex) / ( pnull(i)/pindex_null) - 1] ;


*
*   ---    Set substitution elasticity
*
 sigma = 4;
 rho   = (sigma-1)/sigma;


*
*   ---    Benchmark price/quantity framework
*
 p(i)    = pnull(i);
 x.fx(i) = xnull(i);



*
*   ---    Price index scaled to average import price
*
 pindex_null = sum(i, pnull(i) * xnull(i)) / sum(i, xnull(i));
 pindex.l    = pindex_null;
 pindex.lo   = eps;



*
*   ---     Utility scaled to total imports (quantity)
*
 u.fx        = sum(i, xnull(i));


*
*                  CES form with Hicks-neutral and augmented technological change
*                  Here we set the technological change parameters to unity
*                  i.e. no technological change is assumed so far
*
 model ces_tchange "CES demand system with augmented technological change parameter" / index_aug, demand_aug/;
*

*
*   ---   The Hicks-neutral preference-change parameter is fixed throughout the calibration/simulation
*
  Hicks.fx      =  1;


*
*   ---   The Hicks non-neutral preference-change parameters are fixed to unity in Step 1
*
  theta.fx(i)   =  1;


*
*   ---    Fixed utility assumption throughout (no income effects are taken into account)
*
  u.fx          =  u.l;

*
*   ---     Initialize variables for the calibration of the CES share parameters
*           Step 1 of the two-step calib. procedure for the KVT approach
*
  delta.l (i)   = 1;
  delta.lo(i)   = eps;
  delta.up(i)   = +inf;

*
*   ---     We fix the scale parameter (gamma) to unity 
*
  gamma.fx      = 1;

 ces_tchange.solprint   = 1;
 solve ces_tchange using CNS;
 if(ces_tchange.numinfes ne 0, abort "problem with calibrating the CES demand system -- infeasibilities");
 if(ces_tchange.numredef ne 0, abort "problem with calibrating the CES demand system -- redefs");

   calib_results("gamma", " ", "tchange","observed") = gamma.l;
   calib_results("delta", i  , "tchange","observed") = delta.l(i);
   calib_results("utility"," ","tchange","observed") = u.l;
   calib_results("pindex", " ","tchange","observed") = pindex.l;
   calib_results("theta", i  ,"tchange","observed")  = theta.l(i);
   calib_results("Hicks", " ", "tchange","observed") = Hicks.l;

  option calib_results:3:2:2;


*
*   --- test if the calibrated modle recovers utility
*
 model test_cesaug /ces_aug/;

  u.lo        = 0;
  u.up        = +inf;
  x.fx(i)     = x.l(i);
  delta.fx(i) = delta.l(i);

  test_cesaug.solprint   = 1;
  test_cesaug.iterlim    = 0;
 solve test_cesaug using CNS;
 if(test_cesaug.suminfes gt 1.E-4, abort "problem with recovering utility -- infeasibility");
 if(test_cesaug.numredef gt 1.E-4, abort "problem with recovering utility -- redefs");



*
*  --- test scenario to get 'standard CES' reactions
*      Later we will compare the reaction with the commitment CES to this benchmark
*
*  - relative price of good decreases significantly
 p("2")       = p("2") * .5;


*
*   ---  free up variables, set appropriate bounds
*
 x.l(i)       = 1;
 x.lo(i)      = eps;
 x.up(i)      = +inf;
 pindex.lo    = eps;
 pindex.up    = +inf;


*
*   ---   fixed utility assumption
*
 u.fx         = u.l;

 parameter results(i,*,*)  "simulated results, reporting parameter";

*
*   ---   Test model
*
 model test_tchange / index_aug, demand_aug, Hicks_neutral/;

*
*   ---   Note that as we fixed the Hicks-neutral techn. change parameter to unity, the gamma param. should not change
*         We check this after solving the test model
*
 gamma.lo     = 0.1 * gamma.l;
 gamma.up     = gamma.l * 10;


 test_tchange.solprint   = 1;
 solve test_tchange using CNS;
 if(test_tchange.numinfes ne 0, abort "problem with test simulation -- infeasibilities");
 if(test_tchange.numredef ne 0, abort "problem with test simulation -- redefs");
 if(abs(gamma.l - 10 * gamma.lo) gt 1.E-3, abort "gamma parameter changed!");

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
*  --- step 2 of the KVT calibration approach: Calibrate the 'thetas' to the expected point
*


* Construct an expected point
* ---------------------------
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

* -- alternatively, simply use a guesstimate for the decrease in import demand
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
*   ---  Note that the calibration is not independent from the first calibration to the observed point!


*
*   ---  Price index is calculated automatically by index_aug
*        Here we only initialize the variable with the new average import price
*
 pindex.l    = sum(i, p(i) * x.l(i)) / sum(i, x.l(i));
 pindex.lo   = eps;
 pindex.up   = 3.0 * pindex.l;

*        -  we keep utility at the prevoius level (stay on the same isoquant)
 u.fx          = u.l;


*
*   ---   We solve for thetas; all other CES parameters are fixed
*
 gamma.fx      = gamma.l;
 delta.fx(i)   = delta.l(i);
 Hicks.fx      =  1;

*
*   ---   Initialize the calibrated variable
*
  theta.l (i)   =  1;
  theta.lo(i)   = eps;
  theta.up(i)   = 10.0 * theta.l(i);

 ces_tchange.solprint  = 1;
 solve ces_tchange using CNS;
 if(ces_tchange.suminfes gt 1.E-4, abort "thetas could not be calibrated to the expected point -- infeasibilities");
 if(ces_tchange.numredef gt 1.E-4, abort "thetas could not be calibrated to the expected point -- redefs");

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
 if(isoquant_modKVT.suminfes gt 1.E-4, abort "problem with isoquant test -- infeasibility");
 if(isoquant_modKVT.numredef gt 1.E-4, abort "problem with isoquant test -- redefs");




*
*   ---  Test simulation with the modified KVT approach
*        a) does it reproduce the observed point?

model mod_KVT /index_aug, demand_aug, transition/;


*
*   --- theta will be defined by the transition function
*
 theta.l (i)   = calib_results("theta", i,   "tchange","observed");
 theta.lo(i)   = eps;
 theta.up(i)   = 10 * theta.l(i);

 gamma.fx      = calib_results("gamma", " ", "tchange","observed");

 x.l (i)       = xnull(i);
 x.lo(i)       = eps;
 x.up(i)       = 10 * x.l(i);

 u.fx          = u.l;


*
*   ---     initial price/quantity framework
*
 p(i)          = pnull(i);
 pindex.l      = pindex_null;
 pindex.lo     = eps;
 pindex.up     = 10 * pindex.l;

 mod_KVT.solprint   = 1;
 mod_KVT.iterlim    = 0;
 solve mod_KVT using CNS;
 if(mod_KVT.numinfes ne 0, abort "problem with replicating benchmark with the MODIFIED KVT -- infeasibility");
 if(mod_KVT.numredef ne 0, abort "problem with replicating benchmark with the MODIFIED KVT -- redef");

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
 if(mod_KVT.numinfes ne 0, abort "problem with replicating benchmark with the MODIFIED KVT -- infeasibility");
 if(mod_KVT.numredef ne 0, abort "problem with replicating benchmark with the MODIFIED KVT -- redef");


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
parameter     p_results(scen, i, *)   "reporting parameter";

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



* Use Gnuplot to generate picture
execute 'call %gnuplot_path%gnuplot plot.plt';

* Use mspaint(Windows) to open image file
execute 'mspaint plot.png';


