#!bin/bash

g++ -g -std=c++11 -o opt main_opt.cpp opt.cpp opt.h

num=500
./opt $num

#gnuplot -c exename filename dataname x1 x2 y1 y2
gnuplot -c test.gp opt_SCH.png SCH 0 4 0 4 SCH_opt.txt

gnuplot -c test.gp opt_FON.png FON 0 1 0 1 FON_opt.txt

# gnuplot -c test.gp opt_KUR.png KUR -20 -14 -12 2 KUR_opt.txt

gnuplot -c test.gp opt_ZDT1.png ZDT1 0 1 0 1 ZDT1_opt.txt

gnuplot -c test.gp opt_ZDT2.png ZDT2 0 1 0 1 ZDT2_opt.txt

gnuplot -c test.gp opt_ZDT3.png ZDT3 0 1 -1 1 ZDT3_opt.txt

gnuplot -c test.gp opt_ZDT4.png ZDT4 0 1 0 2 ZDT4_opt.txt

gnuplot -c test.gp opt_ZDT6.png ZDT6 0.2 1 0 1 ZDT6_opt.txt

gnuplot -c test.gp opt_UF1.png UF1 0 1.2 0 1.2 UF1_opt.txt

gnuplot -c test.gp opt_UF2.png UF2 0 1.2 0 1.2 UF2_opt.txt

gnuplot -c test.gp opt_UF3.png UF3 0 1.2 0 1.2 UF3_opt.txt

gnuplot -c test.gp opt_UF4.png UF4 0 1.2 0 1.2 UF4_opt.txt

gnuplot -c test.gp opt_UF5.png UF5 0 1.2 0 1.2 UF5_opt.txt

gnuplot -c test.gp opt_UF6.png UF6 0 1.2 0 1.2 UF6_opt.txt

gnuplot -c test.gp opt_UF7.png UF7 0 1.2 0 1.2 UF7_opt.txt
