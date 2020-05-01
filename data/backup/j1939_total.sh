#!/bin/gnuplot --persist

set log x
set log y
set yrange[1e-17:5e-11]
plot 'j1939_total.dat'
