#!/bin/bash
#awk '{x=$1/100;x=x>0?int(x+0.5):int(x-0.5);++a[x]}END{for (x in a) print 100*x,a[x]}' measurements > measurement_histogram100 
gnuplot << EOF
set terminal postscript eps size 5, 3 enhanced color \
    font 'Helvetica,20' linewidth 2
set output "measurement_histogram.eps"
u(x)=x>1?u(log(x))/(5*x):(1+(1-(2*x-1)**2)/8)/5;
set xlabel "intensity";
set ylabel "frequency";
g=170;
a=1.0/16;
set xrange [-2000:10000];
set logscale y
p "measurement_histogram10" u 1:(\$2/10) w histeps lc rgb "black" title "measurements",\
  418620054/(pi*g*(1+(x/g)**2)) lc rgb "blue" title "outlier distribution"
EOF
