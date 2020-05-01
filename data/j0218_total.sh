#!/bin/gnuplot --persist

set log x
set log y
set yrange[1e-17:5e-11]
plot 'j0218_total.dat' with lines, 'j0218_sync.dat' with lines, 'j0218_IC.dat' with lines,\
    'j0218_cur.dat' with lines
