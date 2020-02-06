#!bin/bash

make clean dep all

run=1
iter=100

#./exe run iter dataset
#gnuplot -c exename filename dataname x1 x2 y1 y2 datafile

./m_nsga2 $run $iter FON 
gnuplot -c gnuplot/test.gp gnuplot/NSGA_FON.png FON 0 1 0 1 pareto/FON/FON_mnsga2_0.txt opt/FON_opt.txt

./m_nsga2 $run $iter SCH
gnuplot -c gnuplot/test.gp gnuplot/NSGA_SCH.png SCH 0 4 0 4 pareto/SCH/SCH_mnsga2_0.txt opt/SCH_opt.txt

./m_nsga2 $run $iter KUR
gnuplot -c gnuplot/test.gp gnuplot/NSGA_KUR.png KUR -20 -14 -12 2 pareto/KUR/KUR_mnsga2_0.txt opt/KUR_opt.txt

./m_nsga2 $run $iter POL
gnuplot -c gnuplot/test.gp gnuplot/NSGA_POL.png POL 0 6 0 30 pareto/POL/POL_mnsga2_0.txt opt/POL_opt.txt

./m_nsga2 $run $iter ZDT1
gnuplot -c gnuplot/test.gp gnuplot/NSGA_ZDT1.png ZDT1 0 1 0 1 pareto/ZDT1/ZDT1_mnsga2_0.txt opt/ZDT1_opt.txt

./m_nsga2 $run $iter ZDT2
gnuplot -c gnuplot/test.gp gnuplot/NSGA_ZDT2.png ZDT2 0 1 0 1 pareto/ZDT2/ZDT2_mnsga2_0.txt opt/ZDT2_opt.txt

./m_nsga2 $run $iter ZDT3
gnuplot -c gnuplot/test.gp gnuplot/NSGA_ZDT3.png ZDT3 0 1 -1 1 pareto/ZDT3/ZDT3_mnsga2_0.txt opt/ZDT3_opt.txt

