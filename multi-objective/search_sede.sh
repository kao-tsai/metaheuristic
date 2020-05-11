
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

# ./search sede $run 100 100 30 100 SCH

#  > pareto/FON/FON_sede_0.txt
#250 100 8 2 6 4 good choice
# 250 100 4 4 6 2 good choice
# ./search sede $run 5000 100 4 5 10 2 FON 0.2 0.033333
# gnuplot -c gnuplot_sede/test.gp gnuplot_sede/SEDE_FON.png FON 0 1 0 1 pareto/FON/sede/FON_sede_0.txt opt_1000/FON_opt.txt

# # 250 500 0.85 0.4
# ./search sede $run 5000 100 4 5 5 2 SCH 0.2 0.033333
# gnuplot -c gnuplot_sede/test.gp gnuplot_sede/SEDE_SCH.png SCH 0 4 0 4 pareto/SCH/sede/SCH_sede_0.txt opt_1000/SCH_opt.txt

# ./search sede $run 5000 100 4 5 10 2 KUR 0.2 0.033333
# gnuplot -c gnuplot_sede/test.gp gnuplot_sede/SEDE_KUR.png KUR -20 -14 -12 2 pareto/KUR/sede/KUR_sede_0.txt opt_1000/KUR_opt.txt

# ./search sede $run 5000 100 4 5 10 2 POL 0.2 0.033333
# gnuplot -c gnuplot_sede/test.gp gnuplot_sede/SEDE_POL.png POL 0 6 0 30 pareto/POL/sede/POL_sede_0.txt opt_1000/POL_opt.txt

# ./search sede $run 250 100 4 5 5 2 ZDT1 0.2 0.001
# gnuplot -c gnuplot_sede/test.gp gnuplot_sede/SEDE_ZDT1.png ZDT1 0 1 0 1 pareto/ZDT1/sede/ZDT1_sede_0.txt opt_1000/ZDT1_opt.txt

# ./search sede $run 250 100 4 5 5 2 ZDT2 0.2 0.001
# gnuplot -c gnuplot_sede/test.gp gnuplot_sede/SEDE_ZDT2.png ZDT2 0 1 0 1 pareto/ZDT2/sede/ZDT2_sede_0.txt opt_1000/ZDT2_opt.txt

# ./search sede $run 250 100 4 5 5 2 ZDT3 0.2 0.01
# gnuplot -c gnuplot_sede/test.gp gnuplot_sede/SEDE_ZDT3.png ZDT3 0 1 -1 1 pareto/ZDT3/sede/ZDT3_sede_0.txt opt_1000/ZDT3_opt.txt

# ./search sede $run 250 100 4 5 5 2 ZDT4 0.1 0.01
# gnuplot -c gnuplot_sede/test.gp gnuplot_sede/SEDE_ZDT4.png ZDT4 0 1 0 1 pareto/ZDT4/sede/ZDT4_sede_0.txt opt_1000/ZDT4_opt.txt

# ./search sede $run 250 100 4 5 5 2 ZDT6 0.2 0.0005
# gnuplot -c gnuplot_sede/test.gp gnuplot_sede/SEDE_ZDT6.png ZDT6 0.2 1 0 1 pareto/ZDT6/sede/ZDT6_sede_0.txt opt_1000/ZDT6_opt.txt

# ./search sede $run 1500 100 4 5 10 2 UF1 0.2 0.01
# gnuplot -c gnuplot_sede/test.gp gnuplot_sede/SEDE_UF1.png UF1 0 1.2 0 1.2 pareto/UF1/sede/UF1_sede_0.txt opt_1000/UF1_opt.txt

# ./search sede $run 1500 100 4 5 10 2 UF2 0.2 0.01
# gnuplot -c gnuplot_sede/test.gp gnuplot_sede/SEDE_UF2.png UF2 0 1.2 0 1.2 pareto/UF2/sede/UF2_sede_0.txt opt_1000/UF2_opt.txt

# ./search sede $run 1000 100 4 5 10 2 UF3 0.2 0.01
# gnuplot -c gnuplot_sede/test.gp gnuplot_sede/SEDE_UF3.png UF3 0 1.2 0 1.2 pareto/UF3/sede/UF3_sede_0.txt opt_1000/UF3_opt.txt

# ./search sede $run 1000 100 4 5 10 2 UF4 0.2 0.01
# gnuplot -c gnuplot_sede/test.gp gnuplot_sede/SEDE_UF4.png UF4 0 1.2 0 1.2 pareto/UF4/sede/UF4_sede_0.txt opt_1000/UF4_opt.txt

# ./search sede $run 1000 100 4 5 10 2 UF5 0.2 0.01
# gnuplot -c gnuplot_sede/test.gp gnuplot_sede/SEDE_UF5.png UF5 0 1.2 0 1.2 pareto/UF5/sede/UF5_sede_0.txt opt_1000/UF5_opt.txt

# ./search sede $run 1000 100 4 5 10 2 UF6 0.2 0.01
# gnuplot -c gnuplot_sede/test.gp gnuplot_sede/SEDE_UF6.png UF6 0 1.2 0 1.2 pareto/UF6/sede/UF6_sede_0.txt opt_1000/UF6_opt.txt

# ./search sede $run 1000 100 4 5 10 2 UF7 0.2 0.01
# gnuplot -c gnuplot_sede/test.gp gnuplot_sede/SEDE_UF7.png UF7 0 1.2 0 1.2 pareto/UF7/sede/UF7_sede_1.txt opt_1000/UF7_opt.txt

./search sede $run 1500 105 4 5 10 2 UF8 0.2 0.001
gnuplot -c gnuplot_sede/dim3.gp gnuplot_sede/SEDE_UF8.png UF8 0 1.2 0 1.2 0 1.2 pareto/UF8/sede/UF8_sede_0.txt opt_1000/UF8_opt.txt

# ./search sede $run 5000 100 4 5 10 2 UF9 0.2 0.001
# gnuplot -c gnuplot_sede/dim3.gp gnuplot_sede/SEDE_UF9.png UF9 0 1.2 0 1.2 0 1.2 pareto/UF9/sede/UF9_sede_0.txt opt_1000/UF9_opt.txt

# ./search sede $run 1500 100 4 5 10 2 UF10 0.2 0.001
# gnuplot -c gnuplot_sede/dim.gp gnuplot_sede/SEDE_UF10.png UF10 0 1.2 0 1.2 pareto/UF10/sede/UF10_sede_0.txt opt_1000/UF10_opt.txt
