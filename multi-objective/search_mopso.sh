
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

# ./search mopso $run 100 100 30 100 evaluate;

# ./search mopso $run 3000 100 30 100 FON
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_FON.png FON 0 1 0 1 pareto/FON/mopso/FON_mopso_0.txt opt/FON_opt.txt

# ./search mopso $run 3000 100 30 100 SCH
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_SCH.png SCH 0 4 0 4 pareto/SCH/mopso/SCH_mopso_0.txt opt/SCH_opt.txt

# ./search mopso $run 3000 100 30 100 KUR
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_KUR.png KUR -20 -14 -12 2 pareto/KUR/mopso/KUR_mopso_0.txt opt/KUR_opt.txt

# ./search mopso $run 3000 100 30 100 POL
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_POL.png POL 0 6 0 30 pareto/POL/mopso/POL_mopso_0.txt opt/POL_opt.txt

#  ./search mopso $run 3000 100 30 100 ZDT1
#  gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_ZDT1.png ZDT1 0 1 0 1 pareto/ZDT1/mopso/ZDT1_mopso_0.txt opt/ZDT1_opt.txt

 ./search mopso $run 3000 100 30 100 ZDT2
gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_ZDT2.png ZDT2 0 1 0 1.5 pareto/ZDT2/mopso/ZDT2_mopso_0.txt opt/ZDT2_opt.txt

# ./search mopso $run 3000 100 30 100 ZDT3
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_ZDT3.png ZDT3 0 1 -1 1 pareto/ZDT3/mopso/ZDT3_mopso_0.txt opt/ZDT3_opt.txt

# ./search mopso $run 3000 100 30 100 ZDT4
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_ZDT4.png ZDT4 0 1 0 10 pareto/ZDT4/mopso/ZDT4_mopso_0.txt opt/ZDT4_opt.txt

# ./search mopso $run 3000 100 30 100 ZDT6
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_ZDT6.png ZDT6 0.2 1 0 1 pareto/ZDT6/mopso/ZDT6_mopso_0.txt opt/ZDT6_opt.txt

# ./search mopso $run 3000 100 30 100 UF1
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_UF1.png UF1 0 1.2 0 1.2 pareto/UF1/mopso/UF1_mopso_0.txt opt/UF1_opt.txt

# ./search mopso $run 3000 100 10 100 UF2
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_UF2.png UF2 0 1.2 0 1.2 pareto/UF2/mopso/UF2_mopso_0.txt opt/UF2_opt.txt

# ./search mopso $run 3000 100 10 100 UF3
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_UF3.png UF3 0 1.2 0 1.2 pareto/UF3/mopso/UF3_mopso_0.txt opt/UF3_opt.txt

# ./search mopso $run 3000 100 10 100 UF4
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_UF4.png UF4 0 1.2 0 1.2 pareto/UF4/mopso/UF4_mopso_0.txt opt/UF4_opt.txt

# ./search mopso $run 3000 100 10 100 UF5
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_UF5.png UF5 0 1.2 0 1.2 pareto/UF5/mopso/UF5_mopso_0.txt opt/UF5_opt.txt

# ./search mopso $run 3000 100 10 100 UF6
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_UF6.png UF6 0 1.2 0 1.2 pareto/UF6/mopso/UF6_mopso_0.txt opt/UF6_opt.txt

# ./search mopso $run 3000 100 10 100 UF7
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_UF7.png UF7 0 1.2 0 1.2 pareto/UF7/mopso/UF7_mopso_0.txt opt/UF7_opt.txt

# ./search mopso $run 3000 50 49 100 UF8
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_UF8.png UF8 0 1.2 0 1.2 pareto/UF8/mopso/UF8_mopso_0.txt opt/UF8_opt.txt

# ./search mopso $run 3000 50 49 500 UF9
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_UF9.png UF9 0 1.2 0 1.2 pareto/UF9/mopso/UF9_mopso_0.txt opt/UF9_opt.txt

# ./search mopso $run 3000 50 49 500 UF10
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_UF10.png UF10 0 1.2 0 1.2 pareto/UF10/mopso/UF10_mopso_0.txt opt/UF10_opt.txt