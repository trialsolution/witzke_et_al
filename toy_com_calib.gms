$ontext
  test different approaches for calibrating the commitement version
  of the modified Armington import demand system



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


* 
*
*   ---    define benchmark input data
* 
 parameter
          pnull(i)           "prices"
          xnull(i)           "demand"
          pindex_null        "price index"
 ;

* set a random price relation for benchmark
 option seed=12234;
 pnull(i)    = uniform(0.5,1);




*
*   ---    Initial demand for consumed goods
*
*testing calibration behaviour if good 1 is small compared to good2
 xnull("1") = 0.2;
 xnull("2") = 0.1;

 display "check initial price/quantity framework", pnull, xnull;


*
*   ---    model setup; parameters, variables and equations
* 

 parameter
           rho               "substit. param."
           sigma             "subst. elasticity"
           p(i)              "observed prices"

           calib_results(*,*,*)  "calibration results -- reporting"
           results(i,*,*)        "simulated results"

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


* dual approach yielding Hicksian demand curves
          demand1(i)         "Hicksian demand curves, version 1"
          index1             "price index, version 1"

 ;

 ces1 ..   u =e= sum(i, delta(i) * (x(i) - mu(i))**rho)**(1/rho);

 demand1(i) $ [x.range(i) or delta.range(i)] .. x(i) - mu(i) =e= delta(i)**sigma * (pindex/p(i))** sigma * u;

 index1 ..     pindex =e= sum(i, delta(i)**sigma * p(i)**(1-sigma))**(1/(1-sigma));



*
*   ---    Model parameterization
*          Exogenous substitution elasticity 
* 
 sigma = 4;
 rho   = (sigma-1)/sigma;





*
* ---  I)  first calibrate the standard CES demand system (commitment terms fixed to zero) 
*

 mu.fx(i)    = 0;




*
*   ---     CES without gamma term 
*           Dual approach; the model consists of a share equation and a price index equation
* 
*      
model calib1 /demand1, index1/;


*
*   ---    initialization for the standard CES calibraiton
* 
 p(i)    = pnull(i);
 x.fx(i) = xnull(i);



* scale utility to total consumption quantity in benchmark
* => the price index will be the average price
 pindex_null = sum(i, pnull(i) * xnull(i)) / sum(i, xnull(i));
 pindex.l    = pindex_null;


 u.fx        = sum(i, xnull(i));

 display "show initial price index and utility level", pindex.l, u.l;

*  initialize bounds
 pindex.lo   = eps;

 delta.l (i)   = 1;
 delta.lo(i)   = eps;


 calib1.solprint  = 1;
 solve calib1 using CNS;

 if(abs(pindex.l - pindex_null) gt 1.E-4, abort "price index not recovered");


*  reporting
  calib_results("gamma", " ", "version1")    = 0;
  calib_results("delta", i  , "version1")    = delta.l(i);
  calib_results("utility"," ","version1")    = u.l;
  calib_results("pindex", " ","version1")    = pindex.l;

* abort "check 1st calibration approach", calib_results;




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



  display "check calibrated parameters", calib_results;

*
*         Using a modified CES function form with commitment-terms
*         The approach allows to deal with the small share problem
*         (and with the problem of zero benchmark demand)
*         Demand reactions of initially small demand shares are therefore magnified
*         compared to the standard CES form approach
*
*    -  Note that the above ces1 equations already contain a commitment term "mu"
*    -  Calibrating the commitment term requires an additional degree of freedom. This will
*       be provided by the extra information on expected demand under different relative prices




*
*   ---     II) test scenario to get a 'standard CES' reaction
*           (Later we will compare the reaction with the commitment CES to this benchmark)
* 
* 
*
*  - relative price of good decreases significantly

* price of good 2 drops to half
 p("2")       = p("2") * .5;

* free quantity variables to get standard CES reaction
 x.lo(i)      = eps;
 x.up(i)      = +inf;


 pindex.lo    = eps;
 pindex.up    = +inf;

* fix utility approach
 u.fx         = u.l;


 delta.fx(i)  = calib_results("delta", i,   "version1");



*
*   ---    same model as above for the calibration, but now the calibrated parameters are fix
* 

 calib1.solprint   = 1;
 solve calib1 using CNS;
 if(calib1.numinfes ne 0, abort "problem with simulation CES version 1");

 results(i, "x", "version1") = x.l(i);




*
*   ---   III.)  Commitment version, see Witzke et al. (2005)
*                Set up modified Armington model, calibrate and solve
*

 set points "calibration points"  / observed, expected/;

 positive variables
          pcom(i,points)          "prices -- commitment version"
          xcom(i,points)          "import demand -- commitment version"
          pindexcom(points)       "price indexes -- commitment version"
          ucom(points)            "utility -- commitment versions"
 ;

 parameter          tangent(points)        "tangent to isoquant -- relative prices for reporting";



 equation
         demand_com(i, points)     "import demand equation"
         index_com(points)         "price index"
         ces_com(points)            "utility aggregator"
 ;

 ces_com(points) $ ucom.range(points) ..

     ucom(points) =e= sum(i, delta(i) * (xcom(i,points) - mu(i))**rho)**(1/rho);

 demand_com(i,points) $ [xcom.range(i,points) or delta.range(i)] ..
* -- to enable working with fixed deltas as well
* demand_com(i,points) $ [xcom.range(i,points) or delta.range(i) or mu.range(i)] ..

     xcom(i,points) - mu(i) =e= delta(i)**sigma * (pindexcom(points)/pcom(i,points))** sigma * ucom(points);

 index_com(points) $ pindexcom.range(points)..

     pindexcom(points) =e= sum(i, delta(i)**sigma * pcom(i,points)**(1-sigma))**(1/(1-sigma));


 model calib_commit "calibration model for the CES commitment version" /index_com, demand_com/;



*
*   ---    Initialization
* 

* we include non-zero commitment term for one commodity only
  mu.fx("1")    = 0;
  mu.lo("2")    = 0;
  mu.up("2")    = +inf;
  mu.l ("2")    = 0.5 * xnull("2");


*  relative prices
 pcom.fx(i, "expected")         = p(i);
 pcom.fx(i, "observed")         = pnull(i);
 tangent(points)               = (-1) * pcom.l("1", points) / pcom.l("2", points);



* assume significantly more demand at the same (lowered) relative price for good 2
* expected demand is initialized with standard CES reaction to price changes
* then increased for good 2
 xcom.fx(i, "expected")         = x.l(i);
 xcom.fx("2", "expected")       = x.l("2");
 xcom.fx("1", "expected")       = x.l("1");

 xcom.fx(i, "observed")         = xnull(i);

* free the quantity variables for good 1 (the one with zero commitment term)
* demand adjusts during calibration
 xcom.lo(i, "expected") $ (not mu.range(i))   = eps;
 xcom.up(i, "expected") $ (not mu.range(i))   = +inf;

 option pcom:3:1:1;
 option xcom:3:1:1;
 display "price/quantity setting", pcom.l, xcom.l;






*
*   ---  share parameters are calibrated simultaneously 
*
 delta.lo(i)  = eps;
 delta.up(i)  = +inf;
 delta.l (i)  = calib_results("delta", i,   "version1") * (p(i)/pnull(i));



*
*   ---    price index will adjust
* 
 pindexcom.lo (points)     = eps;
 pindexcom.l  (points)     = calib_results("pindex", " ",   "version1") ;


*  still fix utility assumption
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




$exit
*
*   --- Prepare plot
*

*  requires a gnuplot local installation 
$setlocal gnuplot_path 'C:\Users\Dev\gnuplot\bin\'



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

set axis   /x,y/;
parameter plotrange(axis);

  plotrange("x") =  xcom.l("2","expected") * 3.;
  plotrange("y") =  xcom.l("1","observed") * 2.;


*
*   --- Prepare GNUPLOT script
*
file pltfile /plot.plt/;
put pltfile;
pltfile.nd=5;
putclose
   'set xlabel "import demand from country 1"'/
   'set ylabel "import demand from country 2"'/
   'set title  "Price SA with the commitment version"'/
*   'set key off'/
   'set xrange [0:', plotrange("x"), ']'/
   'set yrange [0:', plotrange("y"), ']'/
   'set term png font arial 13'/
   'set output "plot.png"'/

* -- add observed and expected points
    'set label at ', xnull("1"), ', ', xnull("2") ,'  "" point pointtype 7 pointsize 1'/
    'set label at ', xcom.l("1","expected"), ', ', xcom.l("2","expected") ,'  "" point pointtype 7 pointsize 1'/

* -- plot the indifference curve
   'plot "plot_commitment.dat" using 1:2 title "indifference curve" with lines, \' /

*  -- add the gradients (showing optimal relative prices) both at the observed and expected points
    tangent("observed"),  '*(x-', xnull("1"), ')  +  ', xnull("2") ,',\ ' /
    tangent("expected"),  '*(x-', xcom.l("1","expected"), ')  +  ', xcom.l("2","expected")

;



* Use Gnuplot to generate picture
execute 'call %gnuplot_path%gnuplot plot.plt';

* Use mspaint(Windows) to open image file
execute 'mspaint plot.png';




