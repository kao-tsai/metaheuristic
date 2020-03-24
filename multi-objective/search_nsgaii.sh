
#!bin/bash
clear
make clean
make

run=10
iter=100
cross_p=0.8
mut_p=0.2
pop_num=500
#./exe run iter dataset
#gnuplot -c exename filename dataname x1 x2 y1 y2 datafile

# ./search bnsga2 $run 100 100 30 100 SCH
# > pareto/FON/FON_bnsga2_0.txt

./search bnsga2 $run 500 100 0.9 0.033333 FON
gnuplot -c gnuplot/test.gp gnuplot/NSGA_FON.png FON 0 1 0 1 pareto/FON/bnsga2/FON_bnsga2_0.txt opt/FON_opt.txt

# 250 500 0.85 0.4
./search bnsga2 $run 500 100 0.9 0.03 SCH
gnuplot -c gnuplot/test.gp gnuplot/NSGA_SCH.png SCH 0 4 0 4 pareto/SCH/bnsga2/SCH_bnsga2_0.txt opt/SCH_opt.txt

# ./search bnsga2 $run 500 100 0.9 0.033333 KUR
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_KUR.png KUR -20 -14 -12 2 pareto/KUR/bnsga2/KUR_bnsga2_0.txt opt/KUR_opt.txt

# ./search bnsga2 $run 500 100 0.9 0.033333 POL
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_POL.png POL 0 6 0 30 pareto/POL/bnsga2/POL_bnsga2_0.txt opt/POL_opt.txt

#  ./search bnsga2 $run 500 100 0.9 0.001 ZDT1
#  gnuplot -c gnuplot/test.gp gnuplot/NSGA_ZDT1.png ZDT1 0 1 0 1 pareto/ZDT1/bnsga2/ZDT1_bnsga2_0.txt opt/ZDT1_opt.txt

#  ./search bnsga2 $run 500 100 0.9 0.001 ZDT2
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_ZDT2.png ZDT2 0 1 0 1 pareto/ZDT2/bnsga2/ZDT2_bnsga2_0.txt opt/ZDT2_opt.txt

# ./search bnsga2 $run 500 100 0.9 0.001 ZDT3
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_ZDT3.png ZDT3 0 1 -1 1 pareto/ZDT3/bnsga2/ZDT3_bnsga2_0.txt opt/ZDT3_opt.txt

# ./search bnsga2 $run 500 100 0.9 0.01 ZDT4
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_ZDT4.png ZDT4 0 1 0 5 pareto/ZDT4/bnsga2/ZDT4_bnsga2_0.txt opt/ZDT4_opt.txt

# ./search bnsga2 $run 500 100 0.9 0.001 ZDT6
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_ZDT6.png ZDT6 0.2 1 0 1 pareto/ZDT6/bnsga2/ZDT6_bnsga2_0.txt opt/ZDT6_opt.txt

./search bnsga2 $run 500 100 0.9 0.01 UF1
gnuplot -c gnuplot/test.gp gnuplot/NSGA_UF1.png UF1 0 1.2 0 1.2 pareto/UF1/bnsga2/UF1_bnsga2_0.txt opt/UF1_opt.txt

./search bnsga2 $run 500 100 0.9 0.01 UF2
gnuplot -c gnuplot/test.gp gnuplot/NSGA_UF2.png UF2 0 1.2 0 1.2 pareto/UF2/bnsga2/UF2_bnsga2_0.txt opt/UF2_opt.txt

./search bnsga2 $run 500 100 0.9 0.01 UF3
gnuplot -c gnuplot/test.gp gnuplot/NSGA_UF3.png UF3 0 1.2 0 1.2 pareto/UF3/bnsga2/UF3_bnsga2_0.txt opt/UF3_opt.txt

./search bnsga2 $run 500 100 0.9 0.01 UF4
gnuplot -c gnuplot/test.gp gnuplot/NSGA_UF4.png UF4 0 1.2 0 1.2 pareto/UF4/bnsga2/UF4_bnsga2_0.txt opt/UF4_opt.txt

./search bnsga2 $run 500 100 0.9 0.01 UF5
gnuplot -c gnuplot/test.gp gnuplot/NSGA_UF5.png UF5 0 1.2 0 1.2 pareto/UF5/bnsga2/UF5_bnsga2_0.txt opt/UF5_opt.txt

./search bnsga2 $run 500 100 0.9 0.01 UF6
gnuplot -c gnuplot/test.gp gnuplot/NSGA_UF6.png UF6 0 1.2 0 1.2 pareto/UF6/bnsga2/UF6_bnsga2_0.txt opt/UF6_opt.txt

./search bnsga2 $run 500 100 0.9 0.01 UF7
gnuplot -c gnuplot/test.gp gnuplot/NSGA_UF7.png UF7 0 1.2 0 1.2 pareto/UF7/bnsga2/UF7_bnsga2_0.txt opt/UF7_opt.txt

# ./search bnsga2 $run 500 100 0.9 0.01 UF8
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_UF8.png UF8 0 1.2 0 1.2 pareto/UF8/bnsga2/UF8_bnsga2_0.txt opt/UF8_opt.txt

# ./search bnsga2 $run 500 100 0.9 0.01 UF9
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_UF9.png UF9 0 1.2 0 1.2 pareto/UF9/bnsga2/UF9_bnsga2_0.txt opt/UF9_opt.txt

# ./search bnsga2 $run 500 100 0.9 0.01 UF10
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_UF10.png UF10 0 1.2 0 1.2 pareto/UF10/bnsga2/UF10_bnsga2_0.txt opt/UF10_opt.txt
