#!/bin/bash

gnuplot -persist <<EOF

reset
set title "ACO"
set xlabel "iteration"
set ylabel "obj value"
set xtics 100
set ytics 5
set terminal svg
set key right
set grid
set output 'ACO/fine_tune_alpha.svg'

plot [0:1024][425:500]'ACO/alpha_0.100000_plot_data.txt' using 1:2 with linespoints pi 100 lw 2 title 'alpha=0.1',\
'ACO/alpha_1.000000_plot_data.txt' using 1:2 with linespoints pi 100 lw 2 title 'alpha=1',\
'ACO/alpha_2.000000_plot_data.txt' using 1:2 with linespoints pi 100 lw 2 title 'alpha=2',\
'ACO/alpha_4.000000_plot_data.txt' using 1:2 with linespoints pi 100 lw 2 title 'alpha=4'

reset
set title "ACO"
set xlabel "iteration"
set ylabel "obj value"
set xtics 100
set ytics 5
set terminal svg
set key right
set grid
set output 'ACO/fine_tune_beta.svg'

plot [0:1024][425:500]'ACO/beta_0.100000_plot_data.txt' with linespoints pi 100 lw 2 title 'beta=0.1',\
'ACO/beta_1.000000_plot_data.txt' using 1:2 with linespoints pi 100 lw 2 title 'beta=1',\
'ACO/beta_2.000000_plot_data.txt' using 1:2 with linespoints pi 100 lw 2 title 'beta=2',\
'ACO/beta_4.000000_plot_data.txt' using 1:2 with linespoints pi 100 lw 2 title 'beta=4'


reset
set title "ACO"
set xlabel "iteration"
set ylabel "obj value"
set xtics 100
set ytics 5
set terminal svg
set key right
set grid
set output 'ACO/fine_tune_rho.svg'

plot [0:1024][425:500]'ACO/rho_0.020000_plot_data.txt' using 1:2 with linespoints pi 100 lw 2 title 'rho=0.02',\
'ACO/rho_0.100000_plot_data.txt' using 1:2 with linespoints pi 100 lw 2 title 'rho=0.1',\
'ACO/rho_0.500000_plot_data.txt' using 1:2 with linespoints pi 100 lw 2 title 'rho=0.5',\
'ACO/rho_0.900000_plot_data.txt' using 1:2 with linespoints pi 100 lw 2 title 'rho=0.9'


EOF
