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

set cbrange [ 0.00000 : 15e-07 ]

set output './F.png'
plot "fimage00.dat" matrix with image

set output './F1.png'
plot "fimage01.dat" matrix with image

set output './F2.png'
plot "fimage02.dat" matrix with image

set cbrange [ 0.00000 :  ]
set output './Q.png'
plot "qimage00.dat" matrix with image

set output './Q1.png'
plot "qimage01.dat" matrix with image

set output './Q2.png'
plot "qimage02.dat" matrix with image

set output './U.png'
plot "uimage00.dat" matrix with image

set output './U1.png'
plot "uimage01.dat" matrix with image

set output './U2.png'
plot "uimage02.dat" matrix with image

