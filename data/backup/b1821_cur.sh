#!/bin/gnuplot --persist

set log x
set log y
set yrange[1e-17:5e-11]
plot '/home/wwc129/Twolayer/data/b1821_cur.dat' with lines;
