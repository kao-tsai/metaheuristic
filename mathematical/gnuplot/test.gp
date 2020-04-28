#!/bin/bash

#gnuplot -persist <<EOF

reset
set title ARG2
set xlabel "Iteration"
set ylabel "obj value"
set xtics 200
set ytics ARG3
set grid
set terminal svg
set key right
set grid
set output ARG1

plot [0:ARG4][ARG5:ARG6]ARG7 using 1:2 with linespoints pt 5 pi 200 lw 2 title 'mfo',\
ARG8 using 1:2 with linespoints pt 3 pi 200 lw 2 title 'pso',\
ARG9 using 1:2 with linespoints pt 1 pi 200 lw 2 title 'ga',\


# EOF
