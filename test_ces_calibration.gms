$ontext
  test different CES functional forms

  1) NOT including productivity multiplier term (gamma)
  2) including productivity multiplier term (gamma)
  3) using a modified CES function form with commitment-terms

  - by including gamma we can force the adding up condition on the deltas
  - in all cases the benchmark utility is scaled to total demand quantity
  - the price index in these cases is scaled to the average price

  References
  ----------

    (standard):  Rutherford, Thomas F. 1995. “CES Demand Functions: Hints and Formulae.”

    (commitment): Witzke, Peter, Marcel Adenäuer, Wolfgang Britz, and Thomas Heckelei. 2005.
                  “Modelling EU Sugar Market Scenarios with a Modified Armington Appoach.” In IATRC Symposium. Sevilla, Spain.
$offtext


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
 xnull("2") = 1.E-2;


 parameter
           rho               "substit. param."
           sigma             "subst. elasticity"
           p(i)              "observed prices"
           calib_delta(i,*)  "calibrated deltas -- reporting parameter"
           calib_gamma(*)    "calibrated gamma  -- reporting parameter"
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

 demand1(i) .. x(i) - mu(i) =e= delta(i)**sigma * (pindex/p(i))** sigma * u;

 demand2(i) .. x(i) - mu(i) =e= (1/gamma) * (delta(i)*gamma)**sigma * (pindex/p(i))** sigma * u;

 index1 ..     pindex =e= sum(i, delta(i)**sigma * p(i)**(1-sigma))**(1/(1-sigma));

 index2 ..     pindex =e= (1/gamma) * sum(i, delta(i)**sigma * p(i)**(1-sigma))**(1/(1-sigma));


* first scale the utility to total consumption quantity in benchmark
* then the price index will be the average price
 model calib1 /demand1, index1/;

 rho   = 0.5;
 sigma = 2;
 mu.fx(i)    = 0;

 p(i)    = pnull(i);
 x.fx(i) = xnull(i);

 pindex_null = sum(i, pnull(i) * xnull(i)) / sum(i, xnull(i));
 pindex.l    = pindex_null;
 u.fx        = sum(i, xnull(i));

 display pindex.l, u.l;

*  initialize
 delta.l (i)   = 1;
 delta.lo(i)   = eps;


 calib1.solprint  = 1;
 solve calib1 using CNS;
* automatic calibration test
 if(abs(pindex.l - pindex_null) gt 1.E-4, abort "price index not recovered");

  calib_delta(i, "version1") = delta.l(i);
  calib_gamma("version1")    = 0;

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


* part b: calibrate a CES aggregator with the gamma term
*         and enforcing the adding-up condition for deltas
 model calib2 /demand2, index2, addup_deltas/;


  u.fx          = u.l;
  delta.lo(i)   = eps;
  delta.up(i)   = +inf;
  gamma.l       = 1;

 calib2.solprint   = 1;
 solve calib2 using CNS;

* automatic calibration test
 if(abs(pindex.l - pindex_null) gt 1.E-4, abort "price index not recovered");

   calib_delta(i, "version2") = delta.l(i);
   calib_gamma("version2")    = gamma.l;


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

 delta.fx(i)  = calib_delta(i,"version1");
 gamma.fx     = calib_gamma("version1");

 model ces_sim1 "simulation model for the first CES approach (without gamma term)" /demand1, index1/;

 ces_sim1.solprint   = 1;
 solve ces_sim1 using CNS;
 if(ces_sim1.numinfes ne 0, abort "problem with simulation CES version 1");

 results(i, "x", "version1") = x.l(i);


* set back demand levels to a default
 x.l(i)       = 1;

 delta.fx(i)  = calib_delta(i,"version2");
 gamma.fx     = calib_gamma("version2");

 model ces_sim2 "simulation model for the second CES approach (with gamma term)" /demand2, index2/;

 ces_sim2.solprint   = 1;
 solve ces_sim2 using CNS;
 if(ces_sim2.numinfes ne 0, abort "problem with simulation CES version 2");

 results(i, "x", "version2") = x.l(i);


 display "check simulated results ", results;



 parameter
          pexp(i)        "expected prices"
          xexp(i)        "expected demand"
 ;

* we include non-zero commitment term for one commodity only
  mu.fx("1")  = 0;

* assume significantly more demand at the same (lowered) relative price for good 2
 pexp(i)     = p(i);
 xexp(i)     = x.l(i);
 xexp("2")   = x.l("2") * 10;

 equation
         isoquant       "the observed point should be on the same isoquant"
 ;

 isoquant ..   u =e= sum(i, delta(i) * (xnull(i) - mu(i))**rho)**(1/rho);


 model calib_commit "calibration model for the CES commitment version" /index1, demand1, isoquant/;

 x.fx(i)  = xexp(i);
 p(i)     = pexp(i);

 mu.lo("2")    = -inf;
 mu.up("2")    = +inf;
 mu.l ("2")    = -1;

 delta.lo(i)  = 0;
 delta.up(i)  = +inf;

 calib_commit.solprint   = 1;
 solve calib_commit using CNS;