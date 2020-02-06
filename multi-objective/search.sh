
#!bin/bash
clear
make clean
make

run=1
iter=100
cross_p=0.8
mut_p=0.2
pop_num=500
#./exe run iter dataset
#gnuplot -c exename filename dataname x1 x2 y1 y2 datafile

# ./search bnsga2 $run 100 500 0.9 0.1 FON > pareto/FON/FON_mnsga2_0.txt
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_FON.png FON 0 1 0 1 pareto/FON/FON_mnsga2_0.txt opt/FON_opt.txt

# #250 500 0.85 0.4
# ./search bnsga2 $run 100 500 0.85 0.25 SCH > pareto/SCH/SCH_mnsga2_0.txt
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_SCH.png SCH 0 4 0 4 pareto/SCH/SCH_mnsga2_0.txt opt/SCH_opt.txt

# ./search bnsga2 $run 100 500 0.9 0.05 KUR > pareto/KUR/KUR_mnsga2_0.txt
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_KUR.png KUR -20 -14 -12 2 pareto/KUR/KUR_mnsga2_0.txt opt/KUR_opt.txt

# ./search bnsga2 $run 100 500 0.8 0.25 POL > pareto/POL/POL_mnsga2_0.txt
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_POL.png POL 0 6 0 30 pareto/POL/POL_mnsga2_0.txt opt/POL_opt.txt

#  ./search bnsga2 $run 250 500 0.85 0.001 ZDT1 > pareto/ZDT1/ZDT1_mnsga2_0.txt
#  gnuplot -c gnuplot/test.gp gnuplot/NSGA_ZDT1.png ZDT1 0 1 0 1 pareto/ZDT1/ZDT1_mnsga2_0.txt opt/ZDT1_opt.txt

#  ./search bnsga2 $run 250 500 0.9 0.001 ZDT2 > pareto/ZDT2/ZDT2_mnsga2_0.txt
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_ZDT2.png ZDT2 0 1 0 1 pareto/ZDT2/ZDT2_mnsga2_0.txt opt/ZDT2_opt.txt

# ./search bnsga2 $run 250 500 0.9 0.001 ZDT3 > pareto/ZDT3/ZDT3_mnsga2_0.txt
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_ZDT3.png ZDT3 0 1 -1 1 pareto/ZDT3/ZDT3_mnsga2_0.txt opt/ZDT3_opt.txt

# ./search bnsga2 $run 1000 100 0.99 0.01 ZDT4 > pareto/ZDT4/ZDT4_mnsga2_0.txt
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_ZDT4.png ZDT4 0 1 0 10 pareto/ZDT4/ZDT4_mnsga2_0.txt opt/ZDT4_opt.txt

# ./search bnsga2 $run 500 500 0.9 0.001 ZDT6 > pareto/ZDT6/ZDT6_mnsga2_0.txt
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_ZDT6.png ZDT6 0.2 1 0 1 pareto/ZDT6/ZDT6_mnsga2_0.txt opt/ZDT6_opt.txt


