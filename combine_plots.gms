*
*   --- Prepare GNUPLOT script
*
file pltfile /plot.plt/;
put pltfile;
putclose
   'set xlabel "import demand from country 1"'/
   'set ylabel "import demand from country 2"'/
   'set title  "Combine the two isoquants"'/
*   'set key off'/
   'set xrange [0:2]'/
   'set yrange [0:2.5]'/
   'set term png font arial 13'/
   'set output "plot_combined.png"'/

   'plot "plot_commitment.dat" using 1:2 title "commitment version" with lines, \' /
   '"plot_multi.dat" using 1:2 title "kuiper tongeren" with lines, \' /
   '"horizontal_commit.dat" using 1:2 title "x2 commit" with lines, "vertical_commit.dat" using 1:2 title "x1 commit" with lines, \' /
   '"horizontal_multi.dat"  using 1:2 title "x2 multi" with lines, "vertical_multi.dat" using 1:2 title "x1 multi" with lines, \'  /
   '"horizontal_obs.dat" using 1:2 title "x2 obs." with lines, "vertical_obs.dat" using 1:2 title "x1 obs." with lines' 
;


* sets and parameters for gnuplotxyz
$setlocal gnuplot_path 'S:\util\gnuplot\bin\'

* Use Gnuplot to generate picture
execute 'call %gnuplot_path%gnuplot plot.plt';

* Use mspaint(Windows) to open image file
execute 'mspaint plot_combined.png';



