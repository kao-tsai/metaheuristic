reset                                                                           
set xlabel 'f_1'
set ylabel 'f_2'
set zlabel 'f_3'
set title ARG2
set term png enhanced font 'Verdana,10'
set key left bottom
set output ARG1

splot[ARG3:ARG4][ARG5:ARG6][ARG7:ARG8]ARG9 using 1:2:3 with points pt 6 lc rgb "web-blue" title 'pareto',\
#ARG8 using 1:2 with lines lw 2 lc rgb "light-red" title 'opt'
#plot [ARG3:ARG4][ARG5:ARG6]ARG7 using 1:2 pt 7 ps 1 title 'pareto'


#plot [:][:500]'iteration.txt' using 1:2 with points title 'iteration',\
#'byte.txt' using 1:2 with points title 'byte',\
#'binary.txt' using 1:2 with points title 'binary',\
#'recursive.txt' using 1:2 with points title 'recursive',\
#'harley.txt' using 1:2 with points title 'harley'




#set term pngcairo font "SetoFont"
#set xlabel "Try it!"
#set ylabel "Hello World!"
#set title "函數 3^x"
#set xrange [0.1:20]
#set output "test.png"

#plot 3**x lw 2
