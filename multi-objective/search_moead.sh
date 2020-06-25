
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

# ./search bnsga2 $run 100 100 30 100 SCH
# > pareto/FON/FON_bnsga2_0.txt

# ./search bnsga2 $run 300 100 0.9 0.033333 FON
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_FON.png FON 0 1 0 1 pareto/FON/bnsga2/FON_bnsga2_0.txt opt/FON_opt.txt

# 1000 500 0.85 0.4
# ./search bnsga2 $run 1000 100 0.9 0.03 SCH
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_SCH.png SCH 0 4 0 4 pareto/SCH/bnsga2/SCH_bnsga2_0.txt opt/SCH_opt.txt

# ./search bnsga2 $run 1000 100 0.9 0.033333 KUR
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_KUR.png KUR -20 -14 -12 2 pareto/KUR/bnsga2/KUR_bnsga2_0.txt opt/KUR_opt.txt

# ./search bnsga2 $run 1000 100 0.9 0.033333 POL
# gnuplot -c gnuplot/test.gp gnuplot/NSGA_POL.png POL 0 6 0 30 pareto/POL/bnsga2/POL_bnsga2_0.txt opt/POL_opt.txt

#  ./search moead $run 250 99 20 1.0 0.001 ZDT1
#  gnuplot -c gnuplot_moead/test.gp gnuplot_moead/MOEAD_ZDT1.png ZDT1 0 1 0 1 pareto/ZDT1/moead/ZDT1_moead_0.txt opt_1000/ZDT1_opt.txt

#  ./search moead $run 250 99 20 1.0 0.001 ZDT2
#  gnuplot -c gnuplot_moead/test.gp gnuplot_moead/MOEAD_ZDT2.png ZDT2 0 1 0 1 pareto/ZDT2/moead/ZDT2_moead_0.txt opt_1000/ZDT2_opt.txt

#  ./search moead $run 250 99 20 1.0 0.001 ZDT3
#  gnuplot -c gnuplot_moead/test.gp gnuplot_moead/MOEAD_ZDT3.png ZDT3 0 1 -1 1 pareto/ZDT3/moead/ZDT3_moead_0.txt opt_1000/ZDT3_opt.txt

#  ./search moead $run 250 99 20 1.0 0.001 ZDT4
#  gnuplot -c gnuplot_moead/test.gp gnuplot_moead/MOEAD_ZDT4.png ZDT4 0 1 0 1 pareto/ZDT4/moead/ZDT4_moead_15.txt opt_1000/ZDT4_opt.txt

# ./search moead $run 250 99 20 1.0 0.001 ZDT6
# gnuplot -c gnuplot_moead/test.gp gnuplot_moead/MOEAD_ZDT6.png ZDT6 0.2 1 0 1 pareto/ZDT6/moead/ZDT6_moead_0.txt opt_1000/ZDT6_opt.txt

# ./search moead $run 3000 99 20 1.0 0.001 UF1
# gnuplot -c gnuplot_moead/test.gp gnuplot_moead/MOEAD_UF1.png UF1 0 1.2 0 1.2 pareto/UF1/moead/UF1_moead_0.txt opt_1000/UF1_opt.txt

 ./search moead $run 500 23 20 1.0 0.001 DTLZ1
 gnuplot -c gnuplot_moead/dim3.gp gnuplot_moead/MOEAD_DTLZ1.png DTLZ1 0 1 0 1 0 1 pareto/DTLZ1/moead/DTLZ1_moead_0.txt opt_1000/DTLZ1_opt.txt

# ./search moead $run 500 23 20 1.0 0.001 DTLZ2
#  gnuplot -c gnuplot_moead/dim3.gp gnuplot_moead/MOEAD_DTLZ2.png DTLZ2 0 1 0 1 0 1 pareto/DTLZ2/moead/DTLZ2_moead_0.txt opt_1000/DTLZ2_opt.txt

