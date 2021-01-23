#!/usr/bin/gnuplot -c

set terminal png font arial 14 size 800,600

set border 3
unset title
unset key
unset xtics
unset ytics
set grid

set xrange [0:]
set yrange [0:]

set cbrange [ 0.00000 :  ]

set output ARG2
plot ARG1 matrix with image
