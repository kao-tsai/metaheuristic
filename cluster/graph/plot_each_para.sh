#!/bin/bash

gnuplot -persist <<EOF

reset
set title "IRIS"
set xlabel "iteration"
set ylabel "obj value"
set xtics 100
set ytics 100
set terminal svg
set key right
set grid
set output 'iris.svg'

plot [0:1024][0:400]'plot_data/ga_iris.txt' using 1:2 with linespoints pi 100 lw 2 title 'ga',\
'plot_data/goa_iris.txt' using 1:2 with linespoints pi 100 lw 2 title 'goa',\
'plot_data/goaga_iris.txt' using 1:2 with linespoints pi 100 lw 2 title 'goaga'

reset
set title "VOWEL"
set xlabel "iteration"
set ylabel "obj value"
set xtics 1000
set ytics 1000
set terminal svg
set key right
set grid
set output 'vowel.svg'

plot [0:10000][1000:5000]'plot_data/ga_vowel.txt' with linespoints pi 1000 lw 2 title 'ga',\
'plot_data/goa_vowel.txt' using 1:2 with linespoints pi 1000 lw 2 title 'goa',\
'plot_data/goaga_vowel.txt' using 1:2 with linespoints pi 1000 lw 2 title 'goaga'


reset
set title "ABALONE"
set xlabel "iteration"
set ylabel "obj value"
set xtics 1000
set ytics 1000
set terminal svg
set key right
set grid
set output 'abalone.svg'

plot [0:10000][4000:20000]'plot_data/ga_abalone.txt' using 1:2 with linespoints pi 1000 lw 2 title 'ga',\
'plot_data/goa_abalone.txt' using 1:2 with linespoints pi 1000 lw 2 title 'goa',\
'plot_data/goaga_abalone.txt' using 1:2 with linespoints pi 1000 lw 2 title 'goaga'

reset
set title "WINE"
set xlabel "iteration"
set ylabel "obj value"
set xtics 100
set ytics 30000000
set terminal svg
set key right
set grid
set output 'wine.svg'

plot [0:1024][15000000:200000000]'plot_data/ga_wine.txt' using 1:2 with linespoints pi 100 lw 2 title 'ga',\
'plot_data/goa_wine.txt' using 1:2 with linespoints pi 100 lw 2 title 'goa',\
'plot_data/goaga_wine.txt' using 1:2 with linespoints pi 100 lw 2 title 'goaga'


EOF
