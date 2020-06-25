
#!bin/bash
clear
make clean
make

run=30
iter=100
cross_p=0.8
mut_p=0.2
pop_num=500
#./exe run iter dataset
#gnuplot -c exename filename dataname x1 x2 y1 y2 datafile

# ./search agmopso $run 25000 99 20 0.9 0.001 SCH
# gnuplot -c gnuplot_agmopso/test.gp gnuplot_agmopso/AGMOPSO_ZDT1.png ZDT1 0 1 0 1 pareto/ZDT1/agmopso/ZDT1_agmopso_0.txt opt_1000/ZDT1_opt.txt

# ./search bnsga2 $run 300 100 0.9 0.033333 FON
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_FON.png FON 0 1 0 1 pareto/FON/bnsga2/FON_bnsga2_0.txt opt/FON_opt.txt

# 1000 500 0.85 0.4
# ./search bnsga2 $run 1000 100 0.9 0.03 SCH
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_SCH.png SCH 0 4 0 4 pareto/SCH/bnsga2/SCH_bnsga2_0.txt opt/SCH_opt.txt

# ./search bnsga2 $run 1000 100 0.9 0.033333 KUR
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_KUR.png KUR -20 -14 -12 2 pareto/KUR/bnsga2/KUR_bnsga2_0.txt opt/KUR_opt.txt

# ./search bnsga2 $run 1000 100 0.9 0.033333 POL
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_POL.png POL 0 6 0 30 pareto/POL/bnsga2/POL_bnsga2_0.txt opt/POL_opt.txt

#  ./search agmopso $run 25000 99 20 0.9 0.001 ZDT1
#  gnuplot -c gnuplot_agmopso/test.gp gnuplot_agmopso/AGMOPSO_ZDT1.png ZDT1 0 1 0 1 pareto/ZDT1/agmopso/ZDT1_agmopso_0.txt opt_1000/ZDT1_opt.txt

#  ./search agmopso $run 25000 99 20 1.0 0.001 ZDT2
#  gnuplot -c gnuplot_agmopso/test.gp gnuplot_agmopso/AGMOPSO_ZDT2.png ZDT2 0 1 0 1 pareto/ZDT2/agmopso/ZDT2_agmopso_0.txt opt_1000/ZDT2_opt.txt

#  ./search agmopso $run 25000 99 20 1.0 0.001 ZDT3
#  gnuplot -c gnuplot_agmopso/test.gp gnuplot_agmopso/AGMOPSO_ZDT3.png ZDT3 0 1 -1 1 pareto/ZDT3/agmopso/ZDT3_agmopso_0.txt opt_1000/ZDT3_opt.txt

#  ./search agmopso $run 25000 99 20 1.0 0.001 ZDT4
#  gnuplot -c gnuplot_agmopso/test.gp gnuplot_agmopso/AGMOPSO_ZDT4.png ZDT4 0 1 0 1 pareto/ZDT4/agmopso/ZDT4_agmopso_15.txt opt_1000/ZDT4_opt.txt

#  ./search agmopso $run 25000 99 20 1.0 0.001 ZDT6
#  gnuplot -c gnuplot_agmopso/test.gp gnuplot_agmopso/AGMOPSO_ZDT6.png ZDT6 0.2 1 0 1 pareto/ZDT6/agmopso/ZDT6_agmopso_0.txt opt_1000/ZDT6_opt.txt

# ./search agmopso $run 300000 99 20 1.0 0.001 UF1
# gnuplot -c gnuplot_agmopso/test.gp gnuplot_agmopso/AGMOPSO_UF1.png UF1 0 1.2 0 1.2 pareto/UF1/agmopso/UF1_agmopso_0.txt opt_1000/UF1_opt.txt

# ./search agmopso $run 300000 99 20 1.0 0.001 UF2
# gnuplot -c gnuplot_agmopso/test.gp gnuplot_agmopso/AGMOPSO_UF2.png UF2 0 1.2 0 1.2 pareto/UF2/agmopso/UF2_agmopso_0.txt opt_1000/UF2_opt.txt

# ./search agmopso $run 300000 99 20 1.0 0.001 UF3
# gnuplot -c gnuplot_agmopso/test.gp gnuplot_agmopso/AGMOPSO_UF3.png UF3 0 1.2 0 1.2 pareto/UF3/agmopso/UF3_agmopso_0.txt opt_1000/UF3_opt.txt

# ./search agmopso $run 300000 99 20 1.0 0.001 UF4
# gnuplot -c gnuplot_agmopso/test.gp gnuplot_agmopso/AGMOPSO_UF4.png UF4 0 1.2 0 1.2 pareto/UF4/agmopso/UF4_agmopso_0.txt opt_1000/UF4_opt.txt

# ./search agmopso $run 300000 99 20 1.0 0.001 UF5
# gnuplot -c gnuplot_agmopso/test.gp gnuplot_agmopso/AGMOPSO_UF5.png UF5 0 1.2 0 1.2 pareto/UF5/agmopso/UF5_agmopso_0.txt opt_1000/UF5_opt.txt

# ./search agmopso $run 300000 99 20 1.0 0.001 UF6
# gnuplot -c gnuplot_agmopso/test.gp gnuplot_agmopso/AGMOPSO_UF6.png UF6 0 1.2 0 1.2 pareto/UF6/agmopso/UF6_agmopso_0.txt opt_1000/UF6_opt.txt

# ./search agmopso $run 300000 99 20 1.0 0.001 UF7
# gnuplot -c gnuplot_agmopso/test.gp gnuplot_agmopso/AGMOPSO_UF7.png UF7 0 1.2 0 1.2 pareto/UF7/agmopso/UF7_agmopso_0.txt opt_1000/UF7_opt.txt

# ./search agmopso $run 315000 13 20 1.0 0.001 UF8
# gnuplot -c gnuplot_agmopso/dim3.gp gnuplot_agmopso/AGMOPSO_UF8.png UF8 0 1 0 1 0 1 pareto/UF8/agmopso/UF8_agmopso_0.txt opt_1000/UF8_opt.txt

# ./search agmopso $run 315000 13 20 1.0 0.001 UF9
# gnuplot -c gnuplot_agmopso/dim3.gp gnuplot_agmopso/AGMOPSO_UF9.png UF9 0 1 0 1 0 1 pareto/UF9/agmopso/UF9_agmopso_0.txt opt_1000/UF9_opt.txt

# ./search agmopso $run 315000 13 20 1.0 0.001 UF10
# gnuplot -c gnuplot_agmopso/dim3.gp gnuplot_agmopso/AGMOPSO_UF10.png UF10 0 1 0 1 0 1 pareto/UF10/agmopso/UF10_agmopso_0.txt opt_1000/UF10_opt.txt

 ./search agmopso $run 52500 13 20 1.0 0.001 DTLZ1
 gnuplot -c gnuplot_agmopso/dim3.gp gnuplot_agmopso/AGMOPSO_DTLZ1.png DTLZ1 0 1 0 1 0 1 pareto/DTLZ1/agmopso/DTLZ1_agmopso_0.txt opt_1000/DTLZ1_opt.txt

./search agmopso $run 52500 13 20 1.0 0.001 DTLZ2
 gnuplot -c gnuplot_agmopso/dim3.gp gnuplot_agmopso/AGMOPSO_DTLZ2.png DTLZ2 0 1.25 0 1.25 0 1.25 pareto/DTLZ2/agmopso/DTLZ2_agmopso_0.txt opt_1000/DTLZ2_opt.txt

./search agmopso $run 52500 13 20 1.0 0.001 DTLZ3
 gnuplot -c gnuplot_agmopso/dim3.gp gnuplot_agmopso/AGMOPSO_DTLZ3.png DTLZ3 0 1 0 1 0 1 pareto/DTLZ3/agmopso/DTLZ3_agmopso_0.txt opt_1000/DTLZ3_opt.txt

# ./search agmopso $run 52500 13 20 1.0 0.001 DTLZ4
#  gnuplot -c gnuplot_agmopso/dim3.gp gnuplot_agmopso/AGMOPSO_DTLZ4.png DTLZ4 0 1.5 0 1.5 0 1.5 pareto/DTLZ4/agmopso/DTLZ4_agmopso_0.txt opt_1000/DTLZ4_opt.txt

./search agmopso $run 52500 13 20 1.0 0.001 DTLZ5
 gnuplot -c gnuplot_agmopso/dim3.gp gnuplot_agmopso/AGMOPSO_DTLZ5.png DTLZ5 0 1 0 1 0 1 pareto/DTLZ5/agmopso/DTLZ5_agmopso_0.txt opt_1000/DTLZ5_opt.txt

./search agmopso $run 52500 13 20 1.0 0.001 DTLZ6
 gnuplot -c gnuplot_agmopso/dim3.gp gnuplot_agmopso/AGMOPSO_DTLZ6.png DTLZ6 0 1.5 0 1.5 0 1.5 pareto/DTLZ6/agmopso/DTLZ6_agmopso_0.txt opt_1000/DTLZ6_opt.txt

./search agmopso $run 52500 13 20 1.0 0.001 DTLZ7
 gnuplot -c gnuplot_agmopso/dim3.gp gnuplot_agmopso/AGMOPSO_DTLZ7.png DTLZ7 0 1.5 0 1.5 0 6 pareto/DTLZ7/agmopso/DTLZ7_agmopso_0.txt opt_1000/DTLZ7_opt.txt