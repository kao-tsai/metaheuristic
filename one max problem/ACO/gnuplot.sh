#!/bin/bash

gnuplot -persist <<EOF

reset
set title "ACO"
set xlabel "iteration"
set ylabel "obj value"
set xtics 100
set ytics 10
set terminal eps
set key right
set grid
set output 'ACO/random_page_result.eps'

plot [0:1500][430:500]'ACO/plot_data.txt' using 1:2 w l lt 1 lw 5 title 'ACO'


EOF
