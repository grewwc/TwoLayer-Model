#!/bin/gnuplot --persist

set log x
set log y
set yrange[1e-17:5e-11]
set xrange[1e-3:4.8e8]
plot 'b1821_total.dat' with lines
