
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

# ./search mopso $run 3000 100 30 200 FON > pareto/FON/FON_mopso_0.txt
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_FON.png FON 0 1 0 1 pareto/FON/FON_mopso_0.txt opt/FON_opt.txt

# ./search mopso $run 3000 100 30 200 SCH > pareto/SCH/SCH_mopso_0.txt
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_SCH.png SCH 0 4 0 4 pareto/SCH/SCH_mopso_0.txt opt/SCH_opt.txt

# ./search mopso $run 3000 100 30 200 KUR > pareto/KUR/KUR_mopso_0.txt
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_KUR.png KUR -20 -14 -12 2 pareto/KUR/KUR_mopso_0.txt opt/KUR_opt.txt

# ./search mopso $run 3000 100 30 200 POL > pareto/POL/POL_mopso_0.txt
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_POL.png POL 0 6 0 30 pareto/POL/POL_mopso_0.txt opt/POL_opt.txt

#  ./search mopso $run 3000 100 30 200 ZDT1 > pareto/ZDT1/ZDT1_mopso_0.txt
#  gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_ZDT1.png ZDT1 0 1 0 1 pareto/ZDT1/ZDT1_mopso_0.txt opt/ZDT1_opt.txt


#  ./search mopso $run 3000 100 20 100 ZDT2 > pareto/ZDT2/ZDT2_mopso_0.txt
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_ZDT2.png ZDT2 0 1 0 1.5 pareto/ZDT2/ZDT2_mopso_0.txt opt/ZDT2_opt.txt

# ./search mopso $run 3000 100 30 200 ZDT3 > pareto/ZDT3/ZDT3_mopso_0.txt
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_ZDT3.png ZDT3 0 1 -1 1 pareto/ZDT3/ZDT3_mopso_0.txt opt/ZDT3_opt.txt

./search mopso $run 3000 500 10 200 ZDT4 > pareto/ZDT4/ZDT4_mopso_0.txt
gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_ZDT4.png ZDT4 0 1 0 2 pareto/ZDT4/ZDT4_mopso_0.txt opt/ZDT4_opt.txt

# ./search mopso $run 100 100 30 200 ZDT6 > pareto/ZDT6/ZDT6_mopso_0.txt
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_ZDT6.png ZDT6 0.2 1 0 1 pareto/ZDT6/ZDT6_mopso_0.txt opt/ZDT6_opt.txt

# ./search mopso $run 5000 50 49 500 UF1 > pareto/UF1/UF1_mopso_0.txt
# gnuplot -c gnuplot_pso/test.gp gnuplot_pso/mopso_UF1.png UF1 0 1.2 0 1.2 pareto/UF1/UF1_mopso_0.txt opt/UF1_opt.txt