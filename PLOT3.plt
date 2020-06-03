#!/usr/bin/gnuplot

set terminal png font arial 14 size 800,600
set border 3
unset title
unset key
set xtics
set ytics
set grid
set xrange [0:200]
set yrange [0:200]

set cbrange [ 0.00000 : 40e-07 ]

set output './F.png'
plot "fimage00_45_00.dat" matrix with image

set palette defined ( 0 "green", 20 "blue", 25 "black",30 "red",50 "yellow"  )

set cbrange [ -2e-7 : 2e-7 ]
set output './Q.png'
plot "qimage00_45_00.dat" matrix with image

set cbrange [ -2e-07 : 2e-07 ]
set output './U.png'
plot "uimage00_45_00.dat" matrix with image
