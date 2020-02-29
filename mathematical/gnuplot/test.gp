#!/bin/bash

#gnuplot -persist <<EOF

reset
set title ARG2
set xlabel "evaluation"
set ylabel "obj value"
set xtics 25000
set ytics 0.01
set terminal svg
set key right
set grid
set output 'compare.svg'

plot [ARG3:ARG4][ARG5:ARG6]ARG7 using 1:2 with linespoints pt 2 pi 100 lw 2 title 'HHO',\
ARG8 using 1:2 with linespoints pi 100 lw 2 title 'SEHHO',\


# EOF
