
gnuplot -persist <<EOF
reset
set title "PATH"
set xlabel "x axis"
set ylabel "y axis"
set xtics 10
set ytics 10
set terminal svg
set key right
set grid
set output 'ACO/path_figure.svg'

plot [0:80][0:80]'ACO/path_data.txt' using 1:2 with linespoints lw 2


reset
set title "ACO"
set xlabel "iteration"
set ylabel "obj value"
set xtics 100
set ytics 5
set terminal svg
set key right
set grid
set output 'ACO/best_result.svg'

plot [0:1024][425:500]'ACO/best_data.txt' using 1:2 with linespoints pi 100 lw 2 title 'ACO',\

EOF
